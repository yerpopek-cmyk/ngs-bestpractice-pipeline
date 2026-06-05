#!/usr/bin/env bash
# =============================================================================
#  steps/01_prepare.sh — Step 1: Data Preparation and Alignment
# =============================================================================
#  Performs the following:
#    1. FASTQ → uBAM conversion (FastqToSam)
#    2. Reference alignment via BwaSpark / BWA-MEM
#    3. Creation of all working directories
#
#  uBAM CONCEPT:
#    Raw FASTQ files are "fragile": the order of reads in paired files can
#    desynchronize upon careless copying/compression.
#    uBAM (unaligned BAM) solves this: it stores pairs as a single binary
#    object + embeds Read Group metadata:
#      --SM  → sample name
#      --RG  → read group identifier
#      --PL  → platform (illumina, PacBio...)
#      --LB  → library
#    This metadata is critical for BQSR and downstream analysis.
# =============================================================================

set -euo pipefail

# ── Source config and utilities ──────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/config.sh"
source "${ROOT_DIR}/lib/utils.sh"

STEP_ID="01_prepare"
init_log "$STEP_ID"

# ── Check checkpoint ────────────────────────────────────────────────────────
if is_done "$STEP_ID"; then
    log_info "Step ${STEP_ID} is already completed. Skipping (delete .checkpoints/${STEP_ID}.done to repeat)."
    exit 0
fi

log_step "01" "Data Preparation and Reference Alignment"

# =============================================================================
#  1.1 CREATE DIRECTORIES
# =============================================================================
log_substep "Creating working directories"

log_theory "Pipeline Directory Structure" \
    "out/ubam/          → unaligned BAMs (uBAM, after FastqToSam)" \
    "out/bwa/           → aligned BAMs (after BWA/BwaSpark)" \
    "out/marks_dup_bam/ → BAMs with marked duplicates" \
    "out/final_bam/     → final BAMs after BQSR (quality corrections applied)" \
    "out/gvcf/          → GVCF files (one per sample, format for joint calling)" \
    "out/fastqc/        → FastQC HTML reports for each BAM" \
    "out/qualimap/      → detailed QC statistics on coverage and alignment" \
    "control/           → file lists (e.g., input.list for CombineGVCFs)" \
    "tmp/               → temporary files for GATK and Spark" \
    "logs/              → log files for each step"

for dir in \
    "${OUT_UBAM}" \
    "${OUT_BWA}" \
    "${OUT_MARKS_DUP}" \
    "${OUT_MARKS_DUP_BAM}" \
    "${OUT_BQSR}" \
    "${OUT_FINAL_BAM}" \
    "${OUT_GVCF}" \
    "${OUT_FASTQC}" \
    "${OUT_QUALIMAP}" \
    "${OUT_MULTIQC}" \
    "${OUT_ANNOTATION}" \
    "${OUT_PHYLO}" \
    "${OUT_KINSHIP}" \
    "${OUT_HAPLOGROUPS}" \
    "${TMP_DIR}" \
    "${LOG_DIR}" \
    "${CONTROL_DIR}" \
    "${CHECKPOINT_DIR}"; do
    ensure_dir "$dir"
done

log_ok "All directories created"

# =============================================================================
#  1.3 CHECK INPUT DATA
# =============================================================================
log_substep "Checking input files"

# Check for FASTQ files
fastq_count=$(count_files "${FASTQ_DIR}/*.fastq.gz")
if [[ "$fastq_count" -eq 0 ]]; then
    log_error "FASTQ files not found in ${FASTQ_DIR}/"
    log_error "Ensure fastq/*.fastq.gz exist."
    exit 1
fi
log_info "FASTQ files found: ${fastq_count}"

# Check reference
require_file "${GENOME_FA}" "Reference genome"
require_file "${GENOME_FA}.fai" "Reference FAI index (run: samtools faidx genome.fa)"
# Dict file may have extension .dict or .fa.dict / .fasta.dict
if [[ -f "${GENOME_FA%.*}.dict" ]]; then
    DICT_FILE="${GENOME_FA%.*}.dict"
elif [[ -f "${GENOME_FA}.dict" ]]; then
    DICT_FILE="${GENOME_FA}.dict"
else
    DICT_FILE="${GENOME_FA%.*}.dict"
fi
require_file "${DICT_FILE}" "Reference Dict file (run: gatk CreateSequenceDictionary)"
require_file "${CHRY_MT_BED}" "Y/MT regions BED file"
require_file "${KNOWN_SITES}" "Known-sites VCF for BQSR"

log_ok "Input files verified"

# =============================================================================
#  1.4 ACTIVATE CONDA ENVIRONMENT
# =============================================================================
log_substep "Activating conda environment: ${CONDA_ENV_GATK}"

log_theory "Conda Environments in NGS Pipelines" \
    "Different tools often require conflicting Python versions and libraries." \
    "Conda creates isolated environments. In this pipeline:" \
    "  gatk4       → GATK 4.x (Java), requires Java 11+" \
    "  QC_fastq    → FastQC, QualiMap, MultiQC, bcftools" \
    "  LineageTracker → Y-haplogroups" \
    "  yleaf       → Y-lineages from BAM" \
    "  haplogrep   → mtDNA haplogroups" \
    "  nextflow    → mtDNA-server (Nextflow workflow for chrM)"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_GATK}" 2>/dev/null || {
    log_warn "conda activate ${CONDA_ENV_GATK} failed. Continuing with current environment."
}

# =============================================================================
#  1.5 FASTQ → uBAM (FastqToSam)
# =============================================================================
log_substep "Converting paired FASTQ → uBAM (FastqToSam)"

log_theory "Why convert FASTQ to uBAM?" \
    "FASTQ is a text format: sequence + ASCII qualities." \
    "Paired FASTQ (R1.fastq.gz + R2.fastq.gz) are two SEPARATE files." \
    "Problem: during transfer/compression, pairs can desynchronize." \
    "uBAM (unaligned BAM) solves this:" \
    "  • Stores R1 and R2 as pairs in a SINGLE file (FLAG: 0x1 = paired)" \
    "  • Embeds Read Group (@RG) — sample and platform metadata" \
    "  • Binary format — safer and more compact than FASTQ" \
    "  • Ready for parallel processing via GATK Spark modules" \
    "" \
    "Key FastqToSam flags:" \
    "  -F1 / -F2     → input R1 and R2 FASTQ" \
    "  -O            → output uBAM" \
    "  --SM {name}   → Sample Name (SM field in @RG)" \
    "  --RG {id}     → Read Group ID (ID field in @RG)" \
    "  --PL illumina → Sequencing platform (PL field in @RG)" \
    "  --LB lib1     → Library (LB field in @RG, important for MarkDuplicates)"

log_info "Running FastqToSam for each FASTQ sample (background parallel)"

for fq in "${FASTQ_DIR}"/*.fastq.gz; do
    [ -e "$fq" ] || continue
    # Skip _2 files to only process per sample base
    if [[ "$fq" == *_2.fastq.gz || "$fq" == *_2_trimmed.fastq.gz ]]; then continue; fi
    
    SAMPLE=$(basename "$fq" | sed -E 's/(_1|_trimmed)?\.fastq\.gz//')
    
    F1="$fq"
    F2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
    F2_TRIMMED="${FASTQ_DIR}/${SAMPLE}_2_trimmed.fastq.gz"
    
    if [ ! -f "${OUT_UBAM}/${SAMPLE}.bam" ]; then
        if [ -f "$F2" ]; then
            log_cmd "FastqToSam (Paired): ${SAMPLE}"
            gatk --java-options -Djava.io.tmpdir="${TMP_DIR}"/ FastqToSam -F1 "$F1" -F2 "$F2" -O "${OUT_UBAM}/${SAMPLE}.bam" --SM "${SAMPLE}" --RG "${SAMPLE}" --PL illumina --LB lib1 &
        elif [ -f "$F2_TRIMMED" ]; then
            log_cmd "FastqToSam (Paired): ${SAMPLE}"
            gatk --java-options -Djava.io.tmpdir="${TMP_DIR}"/ FastqToSam -F1 "$F1" -F2 "$F2_TRIMMED" -O "${OUT_UBAM}/${SAMPLE}.bam" --SM "${SAMPLE}" --RG "${SAMPLE}" --PL illumina --LB lib1 &
        else
            log_cmd "FastqToSam (Single): ${SAMPLE}"
            gatk --java-options -Djava.io.tmpdir="${TMP_DIR}"/ FastqToSam -F1 "$F1" -O "${OUT_UBAM}/${SAMPLE}.bam" --SM "${SAMPLE}" --RG "${SAMPLE}" --PL illumina --LB lib1 &
        fi
    fi
done
wait

log_ok "FastqToSam completed"

# =============================================================================
#  1.6 REFERENCE ALIGNMENT (BwaSpark / BWA-MEM)
# =============================================================================
log_substep "Aligning uBAM to reference (BWA-MEM)"

log_theory "BWA-MEM and Alignment Algorithm" \
    "BWA-MEM (Burrows-Wheeler Aligner) is the gold standard for alignment." \
    "The Burrows-Wheeler algorithm transforms the reference into a compressed structure (BWT)," \
    "allowing it to find exact matches of a read suffix in O(1) time." \
    "" \
    "BwaSpark — GATK wrapper for BWA-MEM on Apache Spark:" \
    "  • Splits BAM into partitions (--bam-partition-size 4MB)" \
    "  • Processes each partition in parallel across ${SPARK_CORES} cores" \
    "  • --spark-master local[${SPARK_CORES}] = local Spark on ${SPARK_CORES} cores" \
    "  • --conf spark.driver.memory=${SPARK_MEMORY} = RAM for Spark driver" \
    "" \
    "Y-chromosome specificity:" \
    "  PAR1/PAR2 (pseudoautosomal regions) are homologous to the X-chromosome." \
    "  Reads mapping to PAR map ambiguously (MAPQ = 0). For proper" \
    "  Y analysis, a hard-masked reference or PAR masking is required."
CHRS=$(awk '{print "\x27"$1"\x27"}' "${CHRY_MT_BED}" | tr '\n' ' ')

BWA_CMD="gatk --java-options \"-Xmx512m -Djava.io.tmpdir=${TMP_DIR}/\" SamToFastq -I {} -FASTQ /dev/stdout -INTERLEAVE true | bwa mem -p -t 4 -R '@RG\\tID:{/.}\\tSM:{/.}\\tPL:illumina\\tLB:lib1' ${GENOME_FA} /dev/stdin | samtools sort -@ 2 -m 250M -o ${TMP_DIR}/{/.}_full.bam - && samtools index ${TMP_DIR}/{/.}_full.bam && samtools view -b ${TMP_DIR}/{/.}_full.bam ${CHRS} > ${OUT_BWA}/{/.}.bam && samtools index ${OUT_BWA}/{/.}.bam && rm ${TMP_DIR}/{/.}_full.bam*"

log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${BWA_CMD}' ::: ${OUT_UBAM}/*.bam"

${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
    "if [ ! -f ${OUT_BWA}/{/.}.bam ] || [ \$(stat -c%s ${OUT_BWA}/{/.}.bam) -lt 1000000 ]; then \
     gatk --java-options \"-Xmx512m -Djava.io.tmpdir=${TMP_DIR}/\" SamToFastq -I {} -FASTQ /dev/stdout -INTERLEAVE true | \
     bwa mem -p -t 4 -R '@RG\\tID:{/.}\\tSM:{/.}\\tPL:illumina\\tLB:lib1' ${GENOME_FA} /dev/stdin | \
     samtools sort -@ 2 -m 250M -o ${TMP_DIR}/{/.}_full.bam - && \
     samtools index ${TMP_DIR}/{/.}_full.bam && \
     samtools view -b ${TMP_DIR}/{/.}_full.bam ${CHRS} > ${OUT_BWA}/{/.}.bam && \
     samtools index ${OUT_BWA}/{/.}.bam && \
     rm ${TMP_DIR}/{/.}_full.bam*; fi" \
    ::: "${OUT_UBAM}"/*.bam

log_ok "BWA-MEM alignment and chrY/chrM extraction completed"

# =============================================================================
#  1.7 ADDITIONAL HARDWARE ACCELERATION INFO
# =============================================================================
log_hardware_tip

# =============================================================================
#  1.8 COMPLETION
# =============================================================================
mark_done "$STEP_ID"

log_info "Step 01 completed. Next step: 02_variant_calling.sh"

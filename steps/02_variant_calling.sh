#!/usr/bin/env bash
# =============================================================================
#  steps/02_variant_calling.sh — Step 2: BAM Processing and Variant Calling
# =============================================================================
#  Performs the following:
#    1. MarkDuplicatesSpark   → mark PCR duplicates
#    2. BaseRecalibrator      → build sequencer error model
#    3. ApplyBQSR             → apply base quality corrections
#    4. HaplotypeCaller       → variant calling (gVCF mode) — Y and MT
#    5. CombineGVCFs          → combine gVCF from all samples
#    6. GenotypeGVCFs         → joint genotyping → final VCF
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/config.sh"
source "${ROOT_DIR}/lib/utils.sh"

STEP_ID="02_variant_calling"
init_log "$STEP_ID"

if is_done "$STEP_ID"; then
    log_info "Step ${STEP_ID} is already completed. Skipping."
    exit 0
fi

log_step "02" "Variant Calling: BQSR + HaplotypeCaller + Joint Genotyping"

# Activate conda environment for GATK4
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_GATK}" 2>/dev/null || \
    log_warn "Failed to activate ${CONDA_ENV_GATK}"

# =============================================================================
#  2.1 MARK DUPLICATES SPARK
# =============================================================================
log_substep "Marking PCR duplicates (MarkDuplicatesSpark)"

log_theory "PCR duplicates and why to mark them" \
    "During library preparation, DNA molecules are amplified via PCR." \
    "Some molecules are copied MORE than others → PCR duplicates." \
    "Problem: duplicates artificially inflate coverage and create false" \
    "single nucleotide variants (PCR artifacts)." \
    "" \
    "MarkDuplicates MARKS (does not remove!) duplicates with the SAM 0x400 flag." \
    "HaplotypeCaller BY DEFAULT IGNORERS reads with this flag." \
    "" \
    "Duplicate criteria: two reads with the exact same 5'-end start position" \
    "and identical sequence → only one remains the 'primary' read." \
    "This is why in FastqToSam we specified --LB lib1: MarkDuplicates" \
    "compares duplicates WITHIN a single library, not across libraries."

MARKDUP_CMD="gatk --java-options \"-Xmx4g -Djava.io.tmpdir=${TMP_DIR}/\" MarkDuplicates -I {} -O ${OUT_MARKS_DUP_BAM}/{/} -M ${OUT_MARKS_DUP}/{/.}.txt --TMP_DIR ${TMP_DIR} --CREATE_INDEX true"

log_info "GNU Parallel -j ${PARALLEL_JOBS}: MarkDuplicates"
log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${MARKDUP_CMD}' ::: ${OUT_BWA}/*.bam"

${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
    "if [ ! -f ${OUT_MARKS_DUP_BAM}/{/} ]; then \
        gatk --java-options \"-Xmx4g -Djava.io.tmpdir=${TMP_DIR}/\" MarkDuplicates \
            -I {} -O ${OUT_MARKS_DUP_BAM}/{/} -M ${OUT_MARKS_DUP}/{/.}.txt --TMP_DIR ${TMP_DIR} --CREATE_INDEX true; \
    else \
        echo 'Skipping MarkDuplicates for {/}, output already exists'; \
    fi" \
    ::: "${OUT_BWA}"/*.bam

log_ok "MarkDuplicatesSpark completed"

# =============================================================================
#  2.2 BASE QUALITY SCORE RECALIBRATION (BQSR) — Step 1: BaseRecalibrator
# =============================================================================
log_substep "BQSR: Building sequencer error model (BaseRecalibrator)"

log_theory "BQSR — Base Quality Score Recalibration" \
    "Sequencers assign a Phred quality score to every nucleotide:" \
    "  Q30 = P(error) = 0.001 (99.9% accuracy)" \
    "  Q20 = P(error) = 0.01  (99% accuracy)" \
    "" \
    "Problem: these qualities are SYSTEMATICALLY over- or under-estimated" \
    "depending on the context: read cycle, preceding dinucleotide," \
    "sequencer tile, etc." \
    "" \
    "BaseRecalibrator builds an empirical model:" \
    "  1. Takes all reads and finds mismatches against the reference" \
    "  2. Excludes positions matching known-sites (dbSNP)" \
    "     — these mismatches are REAL variants, not errors!" \
    "  3. Remaining mismatches = sequencer errors" \
    "  4. Builds a model: reported_Q → empirical_Q for every context" \
    "" \
    "Key flags:" \
    "  --known-sites → VCF with known SNPs (dbSNP): positions that" \
    "                  are NOT considered sequencer errors"

BQSR_TABLE_CMD="gatk --java-options -Djava.io.tmpdir=${TMP_DIR}/ BaseRecalibrator \
    -I {} \
    -R ${GENOME_FA} \
    --known-sites ${KNOWN_SITES} \
    -O ${OUT_BQSR}/{/.}.table \
    --tmp-dir ${TMP_DIR}"

log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${BQSR_TABLE_CMD}' ::: ${OUT_MARKS_DUP_BAM}/*.bam"

${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
    "if [ ! -f ${OUT_BQSR}/{/.}.table ]; then \
        gatk --java-options -Djava.io.tmpdir=${TMP_DIR}/ BaseRecalibrator \
            -I {} -R ${GENOME_FA} \
            --known-sites ${KNOWN_SITES} \
            -O ${OUT_BQSR}/{/.}.table --tmp-dir ${TMP_DIR}; \
    else \
        echo 'Skipping BaseRecalibrator for {/}, output already exists'; \
    fi" \
    ::: "${OUT_MARKS_DUP_BAM}"/*.bam

log_ok "BaseRecalibrator completed"

# =============================================================================
#  2.3 BQSR — Step 2: ApplyBQSR
# =============================================================================
log_substep "BQSR: Applying quality corrections (ApplyBQSR)"

log_theory "ApplyBQSR — applying the correction model" \
    "ApplyBQSR takes the BAM and the table from BaseRecalibrator," \
    "and rewrites Phred qualities in the BAM according to the model." \
    "Result: final_bam/*.bam — the 'golden' file for variant calling." \
    "" \
    "After BQSR, the Q-score distribution more honestly reflects the actual" \
    "error probability. This is critical for HaplotypeCaller:" \
    "HC's Bayesian statistics directly use Phred qualities" \
    "when computing genotype probabilities."

APPLY_BQSR_CMD="gatk --java-options -Djava.io.tmpdir=${TMP_DIR}/ ApplyBQSR \
    -R ${GENOME_FA} \
    -I {} \
    --bqsr-recal-file ${OUT_BQSR}/{/.}.table \
    -O ${OUT_FINAL_BAM}/{/.}.bam \
    --tmp-dir ${TMP_DIR}"

log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${APPLY_BQSR_CMD}' ::: ${OUT_MARKS_DUP_BAM}/*.bam"

    ${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
        "if [ ! -f ${OUT_FINAL_BAM}/{/.}.bam ]; then \
            gatk --java-options -Djava.io.tmpdir=${TMP_DIR}/ ApplyBQSR \
                -R ${GENOME_FA} -I {} \
                --bqsr-recal-file ${OUT_BQSR}/{/.}.table \
                -O ${OUT_FINAL_BAM}/{/.}.bam --tmp-dir ${TMP_DIR}; \
        else \
            echo 'Skipping ApplyBQSR for {/}, output already exists'; \
        fi" \
        ::: "${OUT_MARKS_DUP_BAM}"/*.bam

log_ok "ApplyBQSR completed. Final BAMs → ${OUT_FINAL_BAM}/"

# =============================================================================
#  2.4 VARIANT CALLING: HaplotypeCaller (gVCF mode)
# =============================================================================
log_substep "Variant calling: HaplotypeCaller (gVCF mode)"

log_theory "HaplotypeCaller: local haplotype assembly algorithm" \
    "HC operates on ACTIVE REGIONS — areas with excess mismatches." \
    "" \
    "Algorithm (simplified):" \
    "  1. ACTIVE ZONE: finds a region with suspicious reads" \
    "  2. DE NOVO ASSEMBLY: builds a de Bruijn graph of all possible haplotypes" \
    "     (accounting for SNPs, indels, and their combinations)" \
    "  3. ALIGNMENT: pairwise alignment of each read against each haplotype" \
    "  4. BAYESIAN GENOTYPING: P(genotype | data) = P(data | genotype) × P(genotype)" \
    "     Uses Phred qualities from BQSR!" \
    "  5. Best genotype → output to VCF" \
    "" \
    "Mode -ERC GVCF (genomic VCF):" \
    "  Records info FOR ALL POSITIONS (both variant and invariant)." \
    "  This allows us to later merge samples and perform joint genotyping." \
    "  Without gVCF, we cannot distinguish 'no variant' from 'no data'." \
    "" \
    "Key flags:" \
    "  -ERC GVCF   → enable gVCF generation mode" \
    "  -L chrY_MT.bed → analyze ONLY Y and MT (saves 99% time)" \
    "  -Xmx16g     → allocate 16 GB Java heap (HC is very RAM hungry)" \
    "" \
    "HaplotypeCaller vs DeepVariant (Google):" \
    "  HC:           Bayes + de Bruijn graph + stats. Great for SNPs/indels." \
    "  DeepVariant:  Convolutional Neural Network (CNN). Converts read pileups into" \
    "                RGB 'images' and classifies genotype as a CV task." \
    "                Better at handling complex structural variations."

HC_CMD='gatk --java-options "-Xmx4g -Djava.io.tmpdir='"${TMP_DIR}"'/" HaplotypeCaller \
    -R '"${GENOME_FA}"' \
    -I {} \
    -O '"${OUT_GVCF}"'/{/.}.gvcf.gz \
    -ERC GVCF \
    --sample-ploidy 1 \
    -L '"${CHRY_MT_BED}"' \
    --tmp-dir '"${TMP_DIR}"

log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${HC_CMD}' ::: ${OUT_FINAL_BAM}/*.bam"

    ${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
        "if [ ! -f ${OUT_GVCF}/{/.}.gvcf.gz ]; then \
            gatk --java-options \"-Xmx4g -Djava.io.tmpdir=${TMP_DIR}/\" HaplotypeCaller \
                -R ${GENOME_FA} -I {} \
                -O ${OUT_GVCF}/{/.}.gvcf.gz \
                -ERC GVCF \
                --sample-ploidy 1 \
                -L ${CHRY_MT_BED} \
                --tmp-dir ${TMP_DIR}; \
        else \
            echo 'Skipping HaplotypeCaller for {/}, output already exists'; \
        fi" \
        ::: "${OUT_FINAL_BAM}"/*.bam

log_ok "HaplotypeCaller completed. gVCFs → ${OUT_GVCF}/"

# =============================================================================
#  2.5 JOINT GENOTYPING: Parallel GenotypeGVCFs + bcftools merge
# =============================================================================
log_substep "Joint Genotyping: Parallel GenotypeGVCFs + bcftools merge"

log_theory "Joint Genotyping Optimization" \
    "Instead of the slow CombineGVCFs step which scales poorly and takes hours," \
    "we run GenotypeGVCFs on each sample individually (in parallel) to produce" \
    "single-sample variant VCF files containing only variant sites." \
    "We then merge these files instantly using bcftools merge." \
    "For haploid Y and MT chromosomes, this yields the same target variant set" \
    "needed for haplogroup/kinship analysis while being 100x faster."

# Step 2.5.1: Create individual VCFs in parallel
log_info "Running GenotypeGVCFs on each sample in parallel..."
mkdir -p "${OUT_DIR}/vcf"

${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
    "if [ ! -f ${OUT_DIR}/vcf/{/.}.vcf.gz ]; then \
        gatk --java-options \"-Xmx4g -Djava.io.tmpdir=${TMP_DIR}/\" GenotypeGVCFs \
            -R ${GENOME_FA} \
            -V {} \
            -O ${OUT_DIR}/vcf/{/.}.vcf.gz \
            -L ${CHRY_MT_BED} \
            --tmp-dir ${TMP_DIR}; \
    else \
        echo 'Skipping GenotypeGVCFs for {/}, output already exists'; \
    fi" \
    ::: "${OUT_GVCF}"/*.gvcf.gz

# Step 2.5.2: Merge individual VCFs into the final multi-sample VCF
log_info "Merging individual VCFs into a single multi-sample VCF using bcftools merge..."

# Ensure indices exist (should be created by GATK, but check to be safe)
for f in "${OUT_DIR}"/vcf/*.vcf.gz; do
    if [ ! -f "${f}.tbi" ] && [ ! -f "${f%.vcf.gz}.vcf.gz.tbi" ]; then
        bcftools index -t "$f"
    fi
done

bcftools merge -O z -o "${OUT_DIR}/chrY_MT_final.vcf.gz" "${OUT_DIR}"/vcf/*.vcf.gz
bcftools index -t "${OUT_DIR}/chrY_MT_final.vcf.gz"

log_ok "Genotyping and Merging completed. Final VCF: ${OUT_DIR}/chrY_MT_final.vcf.gz"

# =============================================================================
#  2.6 BCFTOOLS STATS — quick VCF statistics
# =============================================================================
log_substep "gVCF Statistics (bcftools stats)"

log_theory "VCF Quality Metrics" \
    "After variant calling, it is important to check baseline metrics:" \
    "  Ti/Tv ratio (Transitions/Transversions):" \
    "    Expected values: WGS ~2.0–2.1, WES ~3.0–3.3" \
    "    Too low Ti/Tv (< 1.5) → many artifacts" \
    "  Het/Hom ratio: expected ~1.5–2.0 for diploids" \
    "    For Y and MT: always haploid (1/0, no heterozygotes)" \
    "  Missing rate: proportion of missing genotypes (< 5% = normal)"

BCFTOOLS_CMD='bcftools stats '"${OUT_GVCF}"'/{/.}.gvcf.gz > '"${OUT_GVCF}"'/{/.}.txt'
log_cmd "parallel --progress -j ${PARALLEL_JOBS_LIGHT} '${BCFTOOLS_CMD}' ::: ${OUT_UBAM}/*.bam"

${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS_LIGHT}" \
    "bcftools stats ${OUT_GVCF}/{/.}.gvcf.gz > ${OUT_GVCF}/{/.}.txt" \
    ::: "${OUT_UBAM}"/*.bam

log_ok "Step 02 completed"
mark_done "$STEP_ID"

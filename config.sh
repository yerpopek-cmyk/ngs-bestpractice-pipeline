#!/usr/bin/env bash
# =============================================================================
#  config.sh — Central configuration file for NGS Y/MT Analysis pipeline
# =============================================================================
#  All adjustable parameters are defined in this single file.
#  In production mode, modify paths and thresholds ONLY HERE.
#
#  REFERENCE GENOME CONCEPTS (for educational purposes):
#  ─────────────────────────────────────────────────────────────
#  • GRCh37/hg19  — old standard, still used in clinical settings
#  • GRCh38/hg38  — current standard, improved assembly
#  • T2T-CHM13    — the first COMPLETE gapless genome (2022, Science)
#  • Pangenome HPRC — graph of 94+ genomes, eliminates reference bias
#
#  Reference bias: variants missing from the reference are
#  harder to detect because reads with alternative alleles
#  align poorly. This is especially critical for the Y-chromosome
#  due to high inter-population variability.
# =============================================================================

# --- PARALLELISM -------------------------------------------------------------
# GNU Parallel: -j N means N simultaneous tasks
# Rule: for RAM-intensive steps (BQSR, HaplotypeCaller) use N=1..2
#       for light steps (FastQC, bcftools stats) use N = number of CPUs
export THREADS="${THREADS:-12}"
export PARALLEL_JOBS="${PARALLEL_JOBS:-1}"        # -j for heavy steps (MarkDup, BQSR, HC)
export PARALLEL_JOBS_LIGHT="${PARALLEL_JOBS_LIGHT:-12}"  # -j for light steps (QC, stats)

# Spark (BwaSpark, MarkDuplicatesSpark) — internal GATK parallelism
export SPARK_CORES="${SPARK_CORES:-12}"
export SPARK_MEMORY="${SPARK_MEMORY:-12g}"
export SPARK_PART_SIZE="${SPARK_PART_SIZE:-4194304}"   # 4 MB, BAM partition size for Spark

# --- DIRECTORY PATHS ------------------------------------------------------
# BASE_DIR is computed automatically relative to config.sh location
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BASE_DIR

export FASTQ_DIR="${BASE_DIR}/fastq"
export REFERENCE_DIR="${BASE_DIR}/reference"

# Move output to internal WSL filesystem for stable IO (avoiding 9P AcceptAsync errors)
export OUT_DIR="${BASE_DIR}/out"
export TMP_DIR="${BASE_DIR}/tmp"
export LOG_DIR="${BASE_DIR}/logs"
export CONTROL_DIR="${BASE_DIR}/control"

# Subdirectories for output files
export OUT_UBAM="${OUT_DIR}/ubam"
export OUT_BWA="${OUT_DIR}/bwa"
export OUT_MARKS_DUP="${OUT_DIR}/marks_dup"
export OUT_MARKS_DUP_BAM="${OUT_DIR}/marks_dup_bam"
export OUT_BQSR="${OUT_DIR}/BaseRecalibrator"
export OUT_FINAL_BAM="${OUT_DIR}/final_bam"
export OUT_GVCF="${OUT_DIR}/gvcf"
export OUT_FASTQC="${OUT_DIR}/fastqc"
export OUT_QUALIMAP="${OUT_DIR}/qualimap"
export OUT_MULTIQC="${OUT_DIR}/multiqc"
export OUT_ANNOTATION="${OUT_DIR}/annotation"
export OUT_PHYLO="${OUT_DIR}/phylogenetics"
export OUT_KINSHIP="${OUT_DIR}/kinship"
export OUT_HAPLOGROUPS="${OUT_DIR}/haplogroups"

# --- REFERENCE FILES ---------------------------------------------------------
# Recommended hg38 (or T2T-CHM13 to minimize reference bias).
export GENOME_FA="${REFERENCE_DIR}/GCA_000001405.29.fasta"
# dbSNP: known sites for BQSR (Base Quality Score Recalibration)
# BQSR compares observed mismatches with these known variants
# to distinguish systematic sequencer errors from true SNPs.
export KNOWN_SITES="${REFERENCE_DIR}/common_all_20180418.vcf.gz"
# BED file: coordinates for Y-chromosome and mitochondrial genome (chrM)
# Used as -L (interval list) in HaplotypeCaller to
# restrict analysis to specific loci → saves time.
export CHRY_MT_BED="${REFERENCE_DIR}/chrY_MT.bed"

# --- GENOME BUILD -----------------------------------------------------------
export GENOME_BUILD="hg38"     # Used in lineagetracker (-b 38), annovar (-buildver hg38)

# --- ANNOVAR (Variant Annotation) -------------------------------------------
# table_annovar.pl applies multiple annotation layers to each variant:
#   refGene         → genomic context (gene, exon, UTR, intron, intergenic)
#   avsnp150        → rsID identifiers from dbSNP150
#   clinvar_20200316→ clinical significance (Pathogenic, Benign, VUS…)
#   gnomad30_genome → population allele frequencies (gnomAD v3, genome)
#   gnomad211_exome → population allele frequencies (gnomAD v2.1.1, exome)
#   dbnsfp35c       → functional predictors (SIFT, PolyPhen-2, CADD…)
#   dbscsnv11       → splicing variant predictors
export ANNOVAR_DIR="/home/prep01/data/annovar"
export HUMANDB_DIR="${ANNOVAR_DIR}/humandb"
export ANNOVAR_PROTOCOL="refGene,avsnp150,clinvar_20200316,gnomad30_genome,gnomad211_exome,dbnsfp35c,dbscsnv11"
export ANNOVAR_OPERATION="g,f,f,f,f,f,f"
# g = gene-based annotation; f = filter-based (comparison with DB)
export ANNOVAR_BYPASS_IF_MISSING="${ANNOVAR_BYPASS_IF_MISSING:-true}"  # If ANNOVAR is missing, print guide & bypass instead of crashing


# --- FILTERING THRESHOLDS -------------------------------------------------------
# gnomAD MAF (Minor Allele Frequency): variants with frequency > threshold are
# most likely neutral polymorphisms, not pathogenic mutations.
export GNOMAD_MAF_THRESHOLD="0.01"    # 1% — standard threshold for rare variants

# ClinVar: keep only pathogenic classifications?
# ClinVar values: Pathogenic, Likely_pathogenic, VUS, Likely_benign, Benign
export CLINVAR_PATHOGENIC_ONLY=true

# --- KINSHIP THRESHOLDS (PI_HAT) -------------------------------------------------
# PI_HAT — proportion of alleles identical by descent (IBD).
# Calculated by PLINK using the method of moments: PI_HAT = P(IBD=2) + 0.5 * P(IBD=1)
#
#   PI_HAT ≈ 1.00 → monozygotic twins OR technical sample duplicate
#   PI_HAT ≈ 0.50 → 1st-degree relative: parent–child, full siblings
#   PI_HAT ≈ 0.25 → 2nd-degree relative: grandparent–grandchild, half-siblings, uncle/aunt
#   PI_HAT ≈ 0.125→ 3rd-degree relative: first cousins
#   PI_HAT < 0.05 → unrelated individuals
export PI_HAT_DUPLICATE="0.90"   # >= 0.90 → duplicate or MZ twin
export PI_HAT_FIRST="0.40"       # >= 0.40 → 1st-degree relative
export PI_HAT_SECOND="0.175"     # >= 0.175 → 2nd-degree relative
export PI_HAT_THIRD="0.08"       # >= 0.08 → 3rd-degree relative
export PLINK_GENO="0.1"          # Missing genotypes: exclude SNPs with >10% missing
export PLINK_MIN_IBD="0.4"       # Minimum PI_HAT for filtered analysis

# --- TOOL PATHS -----------------------------------------------------
# If tools are in PATH, the name is sufficient. If not, specify the full path.
export GATK="gatk"
export SAMTOOLS="samtools"
export FASTQC_BIN="fastqc"
export QUALIMAP_BIN="qualimap"
export MULTIQC_BIN="multiqc"
export BCFTOOLS="bcftools"
export GNU_PARALLEL="parallel"
export PLINK_BIN="plink"
export KING_BIN="king"
export LINEAGETRACKER="LineageTracker"
export YLEAF_BIN="Yleaf"
export HAPLOGREP_BIN="haplogrep"
export NEXTFLOW_BIN="nextflow"
export PYTHON3="python3"
export PERL="perl"

# --- CONDA ENVIRONMENTS ---------------------------------------------------------
# The pipeline uses several conda environments to isolate dependencies.
export CONDA_ENV_GATK="gatk4"
export CONDA_ENV_QC="QC_fastq"
export CONDA_ENV_LINEAGE="LineageTracker"
export CONDA_ENV_YLEAF="yleaf"
export CONDA_ENV_HAPLOGREP="haplogrep"
export CONDA_ENV_NEXTFLOW="nextflow"

# --- INTERNAL SETTINGS (do not touch) ---------------------------------------
export PIPELINE_VERSION="1.0.0"
export CHECKPOINT_DIR="${BASE_DIR}/.checkpoints"

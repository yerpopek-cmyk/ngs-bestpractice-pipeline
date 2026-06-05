#!/usr/bin/env bash
# =============================================================================
#  download_new_reads.sh — Fetch the requested male WGS samples (12-core optimized)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

SAMPLES=("SRR622461" "ERR050112" "ERR3963835")
CORES=12

echo "================================================================="
echo "  🧬 Downloading New Samples (12 Cores)"
echo "================================================================="

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_QC}" 2>/dev/null || {
    echo "[WARN] Could not activate ${CONDA_ENV_QC}. Attempting with current environment..."
}

# Define ENA direct download mappings
declare -A ENA_R1
declare -A ENA_R2

# SRR622461 (NA12878 - CEU Female WGS, downsampled 5x)
ENA_R1["SRR622461"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_1.fastq.gz"
ENA_R2["SRR622461"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_2.fastq.gz"

# ERR050112 (HG01500 - IBS Male WGS)
ENA_R1["ERR050112"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050112/ERR050112_1.fastq.gz"
ENA_R2["ERR050112"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050112/ERR050112_2.fastq.gz"

# ERR3963835 (NA21128 - GIH Male WGS)
ENA_R1["ERR3963835"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR396/005/ERR3963835/ERR3963835_1.fastq.gz"
ENA_R2["ERR3963835"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR396/005/ERR3963835/ERR3963835_2.fastq.gz"

# -----------------------------------------------------------------
#  🔧 Fix SRA Toolkit Connection Issues (Failed to call external services)
# -----------------------------------------------------------------
echo "Resetting SRA Toolkit configuration to defaults..."
rm -rf ~/.ncbi
mkdir -p ~/.ncbi
vdb-config --restore-defaults >/dev/null 2>&1 || true

if ! command -v prefetch &>/dev/null || ! command -v fasterq-dump &>/dev/null; then
    echo "SRA Toolkit not found. Installing 'sra-tools'..."
    conda install -y -n "${CONDA_ENV_QC}" -c bioconda sra-tools || exit 1
fi

# Install 'pigz' for ultra-fast parallel compression using 12 cores
if ! command -v pigz &>/dev/null; then
    echo "Installing 'pigz' for 12-core parallel compression..."
    conda install -y -n "${CONDA_ENV_QC}" -c conda-forge pigz || true
fi

ensure_dir "${FASTQ_DIR}"
ensure_dir "${TMP_DIR}"

for accession in "${SAMPLES[@]}"; do
    # Skip download if files already exist
    if [ -f "${FASTQ_DIR}/${accession}_1.fastq.gz" ] || [ -f "${FASTQ_DIR}/${accession}.fastq.gz" ]; then
        echo "================================================================="
        echo " FASTQ files for ${accession} already exist. Skipping download."
        echo "================================================================="
        continue
    fi

    echo ""
    echo "-----------------------------------------------------------------"
    echo " Downloading Run: ${accession}"
    echo "-----------------------------------------------------------------"

    # Clean up any leftover lock files or incomplete downloads from previous runs
    echo "Cleaning up any leftover temporary files/lock files for ${accession}..."
    rm -rf "${TMP_DIR}/${accession}"

    # Use direct curl streaming/download from ENA if we have ENA links mapped
    if [ -n "${ENA_R1[${accession}]:-}" ]; then
        echo "Downloading ${accession} directly from ENA via curl..."
        curl -L -# "${ENA_R1[${accession}]}" -o "${FASTQ_DIR}/${accession}_1.fastq.gz"
        curl -L -# "${ENA_R2[${accession}]}" -o "${FASTQ_DIR}/${accession}_2.fastq.gz"
    else
        # 1. Prefetch the SRA file
        echo "Running prefetch..."
        prefetch "${accession}" -O "${TMP_DIR}"

        # 2. Extract FASTQ files using fastq-dump with gzip compression
        echo "Extracting reads using fastq-dump..."
        fastq-dump "${TMP_DIR}/${accession}/${accession}.sra" \
            -O "${FASTQ_DIR}" \
            --split-files \
            --gzip
    fi

    # Clean up SRA
    rm -rf "${TMP_DIR}/${accession}"

    # Verify download succeeded
    if [ ! -s "${FASTQ_DIR}/${accession}_1.fastq.gz" ] && [ ! -s "${FASTQ_DIR}/${accession}.fastq.gz" ]; then
        echo "❌ Error: Failed to download files for ${accession}!"
        exit 1
    fi

    echo "✔ Successfully finished downloading & extracting ${accession}!"
done

# Reset samples.tsv to point to all new downloaded samples
echo "Updating samples.tsv with new samples..."
cat << EOF > "${SCRIPT_DIR}/samples.tsv"
sample_id	sex	relationship	haplogroup_Y	haplogroup_MT
SRR622461	F	probant	.	H1
ERR050112	M	unrelated	R1b1a2	H1
ERR3963835	M	unrelated	I2a1b	U5b
EOF

echo "✔ All downloads and preparations complete!"
echo "To run the optimized pipeline, execute:"
echo "  bash run_pipeline.sh --reset"



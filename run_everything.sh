#!/usr/bin/env bash
# =============================================================================
#  run_everything.sh — Download raw data, extract, and execute pipeline
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

# Ensure directories exist
ensure_dir "${FASTQ_DIR}"
ensure_dir "${TMP_DIR}"

# 1. Download raw reads using prefetch
echo "================================================================="
echo " 1. Downloading SRA Archives (prefetch)"
echo "================================================================="
# Activate conda environment for QC/SRA toolkit
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_QC}" 2>/dev/null || true

# Remove any leftover lock files from previous runs/crashes
for srr in SRR1976948 SRR1980220 SRR1975755; do
    if [ -f "${TMP_DIR}/${srr}/${srr}.sra.lock" ]; then
        echo "Removing leftover lock file for ${srr}..."
        rm -f "${TMP_DIR}/${srr}/${srr}.sra.lock"
    fi
done

# Parallel downloads to maximize speed
echo "Starting downloads in parallel..."
prefetch SRR1976948 -O "${TMP_DIR}" &
pid1=$!
prefetch SRR1980220 -O "${TMP_DIR}" &
pid2=$!
prefetch SRR1975755 -O "${TMP_DIR}" &
pid3=$!

# Wait for all background downloads to finish
wait $pid1
wait $pid2
wait $pid3
echo "All downloads complete!"

# 2. Extract and Compress FASTQ files
echo "================================================================="
echo " 2. Extracting & compressing FASTQ files (fasterq-dump + pigz)"
echo "================================================================="
for srr in SRR1976948 SRR1980220 SRR1975755; do
    # Skip if compressed fastq files already exist
    if [ -f "${FASTQ_DIR}/${srr}_1.fastq.gz" ] && [ -f "${FASTQ_DIR}/${srr}_2.fastq.gz" ]; then
        echo "✔ FASTQ files for ${srr} already exist. Skipping extraction."
        continue
    fi

    echo "Extracting ${srr}..."
    fasterq-dump --split-files -e "${THREADS}" --mem 8GB -O "${FASTQ_DIR}" "${TMP_DIR}/${srr}/${srr}.sra"

    echo "Compressing ${srr} reads..."
    pigz -f -p "${THREADS}" "${FASTQ_DIR}/${srr}_1.fastq" &
    pigz -f -p "${THREADS}" "${FASTQ_DIR}/${srr}_2.fastq" &
    wait
done

# Cleanup SRA files to save space
echo "Cleaning up temporary SRA files..."
rm -rf "${TMP_DIR}"/SRR*

# 3. Update samples.tsv just in case
echo "Updating samples.tsv..."
cat << EOF > "${SCRIPT_DIR}/samples.tsv"
sample_id	sex	relationship	haplogroup_Y	haplogroup_MT
SRR1976948	M	probant	R1b1a2	H1
SRR1980220	M	unrelated	R1b1a2	H1
SRR1975755	M	unrelated	I2a1b	U5b
EOF

# 4. Run the main pipeline
echo "================================================================="
echo " 3. Running NGS Pipeline (run_pipeline.sh)"
echo "================================================================="
cd "${SCRIPT_DIR}"
# Run pipeline with reset to start clean
bash run_pipeline.sh --reset

echo "================================================================="
echo " 🎉 EVERYTHING COMPLETED SUCCESSFULLY!"
echo "================================================================="

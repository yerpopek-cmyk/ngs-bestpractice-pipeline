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

# Map SRR accession to Sample Name
declare -A SAMPLE_MAP
SAMPLE_MAP["SRR2052353"]="HG001"
SAMPLE_MAP["SRR1766555"]="HG002"
SAMPLE_MAP["SRR1766749"]="HG003"

# Remove any leftover lock files from previous runs/crashes
for srr in SRR2052353 SRR1766555 SRR1766749; do
    if [ -f "${TMP_DIR}/${srr}/${srr}.sra.lock" ]; then
        echo "Removing leftover lock file for ${srr}..."
        rm -f "${TMP_DIR}/${srr}/${srr}.sra.lock"
    fi
done

# Parallel downloads to maximize speed
echo "Starting downloads in parallel..."
prefetch SRR2052353 -O "${TMP_DIR}" &
pid1=$!
prefetch SRR1766555 -O "${TMP_DIR}" &
pid2=$!
prefetch SRR1766749 -O "${TMP_DIR}" &
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
for srr in SRR2052353 SRR1766555 SRR1766749; do
    sample="${SAMPLE_MAP[$srr]}"
    # Skip if compressed fastq files already exist
    if [ -f "${FASTQ_DIR}/${sample}_1.fastq.gz" ] && [ -f "${FASTQ_DIR}/${sample}_2.fastq.gz" ]; then
        echo "✔ FASTQ files for ${sample} already exist. Skipping extraction."
        continue
    fi

    echo "Extracting ${srr} (as ${sample})..."
    fasterq-dump --split-files -e "${THREADS}" --mem 8GB -O "${FASTQ_DIR}" "${TMP_DIR}/${srr}/${srr}.sra"

    # Rename extracted files to the standard HG001/HG002/HG003 names
    mv "${FASTQ_DIR}/${srr}_1.fastq" "${FASTQ_DIR}/${sample}_1.fastq"
    mv "${FASTQ_DIR}/${srr}_2.fastq" "${FASTQ_DIR}/${sample}_2.fastq"

    echo "Compressing ${sample} reads..."
    pigz -f -p "${THREADS}" "${FASTQ_DIR}/${sample}_1.fastq" &
    pigz -f -p "${THREADS}" "${FASTQ_DIR}/${sample}_2.fastq" &
    wait
done

# Cleanup SRA files to save space
echo "Cleaning up temporary SRA files..."
rm -rf "${TMP_DIR}"/SRR*

# 3. Update samples.tsv just in case
echo "Updating samples.tsv..."
cat << EOF > "${SCRIPT_DIR}/samples.tsv"
sample_id	sex	relationship	haplogroup_Y	haplogroup_MT
HG001	F	probant	.	H13a1a1a
HG002	M	son	J1a2a1a2c1a1a1~	H5a7
HG003	M	father	J1a2a1a2c1a1	K1a1b1a
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

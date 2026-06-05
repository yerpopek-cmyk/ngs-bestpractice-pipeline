#!/usr/bin/env bash
# =============================================================================
#  download_quick_test.sh — Fetch a tiny downsampled subset for rapid testing
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

ACCESSION="ERR050112"
NUM_READS=1000000 # 1 million reads (4 million lines in FASTQ)
NUM_LINES=$((NUM_READS * 4))

echo "================================================================="
echo "  🧬 Creating Fast Downsampled Test Dataset (1M reads)"
echo "  Target Sample: ${ACCESSION} (Male, Iberian, Spain)"
echo "  Method: Direct HTTPS Streaming from ENA (EBI)"
echo "================================================================="

ensure_dir "${FASTQ_DIR}"
ensure_dir "${TMP_DIR}"

# Clean up any existing FASTQ/SRA files for this accession in the FASTQ folder to start clean
rm -f "${FASTQ_DIR}/${ACCESSION}"*.fastq*

echo ""
echo "-----------------------------------------------------------------"
echo " Streaming and downsampling reads directly from ENA (Europe)..."
echo " (Bypassing SRA Toolkit completely. Finishes in ~30 seconds!)"
echo "-----------------------------------------------------------------"

# Disable pipefail temporarily to allow head to close the pipe (causing SIGPIPE) without crashing the script
set +o pipefail

echo "Streaming Forward Reads (R1)..."
curl -L -# "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050112/ERR050112_1.fastq.gz" | \
    gunzip -c 2>/dev/null | \
    head -n ${NUM_LINES} | \
    gzip -1 > "${FASTQ_DIR}/${ACCESSION}_1.fastq.gz"

echo "Streaming Reverse Reads (R2)..."
curl -L -# "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050112/ERR050112_2.fastq.gz" | \
    gunzip -c 2>/dev/null | \
    head -n ${NUM_LINES} | \
    gzip -1 > "${FASTQ_DIR}/${ACCESSION}_2.fastq.gz"

set -o pipefail

# Verify both files were created and are larger than 0 bytes
if [ -s "${FASTQ_DIR}/${ACCESSION}_1.fastq.gz" ] && [ -s "${FASTQ_DIR}/${ACCESSION}_2.fastq.gz" ]; then
    echo ""
    echo "✔ Successfully downloaded and downsampled synchronized paired-end reads!"
else
    echo ""
    echo "❌ Error: Streaming failed or files are empty!"
    exit 1
fi

echo "-----------------------------------------------------------------"
echo " Updating samples.tsv for the quick test..."
echo "-----------------------------------------------------------------"
cat << EOF > "${SCRIPT_DIR}/samples.tsv"
sample_id	sex	relationship	haplogroup_Y	haplogroup_MT
${ACCESSION}	M	probant	R1b1a2	H1
EOF

echo "✔ samples.tsv updated with ${ACCESSION}!"
echo ""
echo "================================================================="
echo " 🎉 QUICK TEST SETUP COMPLETE!"
echo " Raw data size: ~100MB instead of 15GB!"
echo " To run the pipeline now, execute:"
echo "   bash run_pipeline.sh --reset"
echo "================================================================="


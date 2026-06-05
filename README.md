# GATK Best Practices Pipeline
######                                           ---NGS Y/MT Analysis Pipeline---

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Pipeline-blue)
![Bash](https://img.shields.io/badge/Language-Bash%20%7C%20Python-green)
![License](https://img.shields.io/badge/License-MIT-orange)

This repository hosts a production-grade bioinformatics pipeline for analyzing both Y-chromosomal (Y-DNA) and Mitochondrial (mtDNA) genomes from Next-Generation Sequencing (NGS) data. The primary focus of this project is reconstructing paternal and maternal lineages, assigning haplogroups, and assessing kinship in the context of the reference genome.

## Experiment Overview

In this experiment, the pipeline was executed on the sequencing sample **SRR622461** against the human reference genome **GCA_000001405.29 (GRCh38)**. 

The pipeline runs through 6 main steps:
1. **Data Preparation & Alignment**: Aligning raw reads using `bwa mem`.
2. **Variant Calling**: Identifying single nucleotide polymorphisms (SNPs) via GATK's `HaplotypeCaller` following base quality recalibration (BQSR).
3. **Quality Control**: Aggregating QC metrics via `FastQC`, `QualiMap`, and `MultiQC`.
4. **Annotation**: Annotating variants for potential pathogenicity (e.g., using ClinVar).
5. **Phylogenetics**: Assigning haplogroups and tracking lineages using specialized tools.
6. **Kinship Analysis**: Determining relatedness between samples.

The output files include comprehensive haplogroup summaries, pathogenic variant annotations, and a MultiQC dashboard, which are included in the `out/` directory of this repository.

## Theoretical Background

Both the Y-chromosome and mitochondrial DNA serve as excellent "molecular clocks" in population genetics and evolutionary biology because they are largely non-recombining. Mutations accumulate linearly over time, allowing us to build precise phylogenetic trees.

### Y-Chromosome (Paternal Lineage)
The Y-chromosome is passed almost exclusively from father to son. Due to its lack of recombination with the X-chromosome (except at the pseudoautosomal regions), it preserves a direct record of paternal ancestry.
* **Lineages**: Different branches of the Y-chromosome tree are termed haplogroups (e.g., R1b, I2a, E1b1a).
* **Tools Used**: `LineageTracker`, `Yleaf`
* **Useful Resources**:
  * [YFull Tree](https://www.yfull.com/tree/)
  * [ISOGG Tree](https://isogg.org/tree/)

### Mitochondrial DNA (Maternal Lineage)
Mitochondrial DNA (mtDNA) is a small circular genome (~16.5 kb) inherited strictly from the mother. It traces maternal lineages back to a common female ancestor.
* **Lineages**: Major haplogroups include H (common in Western Europe), U, and L (ancient African lineages).
* **Tools Used**: `HaploGrep`
* **Useful Resources**:
  * [PhyloTree](https://www.phylotree.org/)
  * [HaploGrep](https://haplogrep.i-med.ac.at/)

## Did You Know? 🧬

* **They Never Met!** The maternal ancestor of all living humans ("Mitochondrial Eve") lived approximately 150,000–200,000 years ago in Africa. The paternal ancestor of all living men ("Y-chromosomal Adam") lived around 200,000–300,000 years ago. Despite their biblical nicknames, they lived tens of thousands of years apart and never met!
* **The Shrinking Y-Chromosome:** Over the last 166 million years, the Y-chromosome has lost most of its original genes, shrinking from over 1,000 genes to just about 55 active protein-coding genes today. Luckily, it seems to have stabilized and isn't going anywhere anytime soon!
* **Ubiquitin Death Sentence:** Why is mitochondrial DNA inherited strictly from the mother? When a sperm fertilizes an egg, the sperm's mitochondria are tagged with a protein called *ubiquitin*, marking them for destruction. The egg's internal cellular machinery promptly destroys them, ensuring only maternal mtDNA persists in the offspring.

## System Requirements & Dependencies

To execute this pipeline, your system must meet the following dependency and reference requirements. All environment builds are handled automatically by the setup script, but you must supply the reference datasets.

### 1. Required Reference Genomes & Databases
* **Reference Genome**: Human Reference Genome **GRCh38 / hg38** (specifically `GCA_000001405.29_GRCh38_no_alt_analysis_set.fna`).
* **Annotation Databases**:
  * **ClinVar**: Used in `04_annotation.sh` for clinical variant significance annotations.
  * **gnomAD**: Used to check population allele frequencies and filter common variants.

### 2. Core Bioinformatic Tools
* **Alignment & Processing**: `bwa` (Burrows-Wheeler Aligner), `samtools`, and `bcftools`.
* **Variant Calling**: `GATK4` (Genome Analysis Toolkit for MarkDuplicates, BQSR, and HaplotypeCaller).
* **QC & Visualization**: `FastQC`, `QualiMap`, and `MultiQC`.
* **Phylogenetics & Lineage Assignment**: `Yleaf` and `LineageTracker` (for Y-DNA haplogroups), `HaploGrep3` (for mtDNA haplogroups).
* **Kinship Analysis**: `PLINK2` and `KING` (for relatedness, coefficient analysis, and pedigree checks).

### 3. Execution & Optimization (GNU Parallel)
This pipeline is engineered for extreme efficiency. It leverages **GNU Parallel** to run jobs concurrently across all available CPU threads:
* **Concurreny Acceleration**: Rather than executing alignment, sorting, duplicate marking, BQSR, and variant calling sequentially, the pipeline splits tasks (across chromosome blocks or sample lanes) and processes them in parallel.
* **Resource Optimization**: CPU thread limits can be fully customized in `config.sh` (via the `PARALLEL_JOBS` and `PARALLEL_JOBS_LIGHT` variables) to perfectly fit your workstation's capacity.

---

## How to Run This Pipeline

To reproduce this experiment, you will need WSL (if on Windows) or a Linux environment with `conda` installed.

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yerpopek-cmyk/ngs-y-mt-pipeline.git
   cd ngs-y-mt-pipeline
   ```

2. **Setup Environments:**
   Run the setup script to initialize the necessary conda environments:
   ```bash
   bash setup_conda_envs.sh
   ```

3. **Provide Raw Data:**
   Place your raw `.fastq` files in the `fastq/` folder, and ensure the reference genome is located in `reference/`.

4. **Execute the Pipeline:**
   Run the complete pipeline in production mode:
   ```bash
   bash run_pipeline.sh --production
   ```
   You can also specify threads or resume from specific steps:
   ```bash
   bash run_pipeline.sh --production --from 02 --threads 8
   ```

## GPU Acceleration with NVIDIA Clara Parabricks

For users with access to high-performance NVIDIA GPUs (requires at least **16 GB VRAM**, e.g., RTX 3090/4090, L4, A100), you can significantly accelerate the alignment and variant-calling steps by replacing the standard CPU-based tools with **NVIDIA Clara Parabricks**. 

This reduces the processing time for alignment and variant calling from hours down to just minutes.

> [!NOTE]
> To run Clara Parabricks, your system must have Docker and the **NVIDIA Container Toolkit** installed and configured to access your GPU.

### 1. GPU-Accelerated Alignment (Replaces BWA-MEM + Samtools + MarkDuplicates)
Instead of running standard BWA-MEM (`steps/01_prepare.sh`), Parabricks' `fq2bam` combines alignment, sorting, duplicate marking, and Base Quality Score Recalibration (BQSR) calculation into a single GPU-optimized step:

```bash
pbrun fq2bam \
  --ref reference/GCA_000001405.29_GRCh38_no_alt_analysis_set.fna \
  --in-fq fastq/SRR622461_1.fastq.gz fastq/SRR622461_2.fastq.gz \
  --out-bam out/alignment/SRR622461_merged.bam \
  --out-recal-file out/alignment/SRR622461_recal.txt \
  --num-gpus 1
```

> [!TIP]
> If you are running on a 16 GB VRAM card (like an RTX 4080 or A2000), add the `--low-memory` flag to ensure the program fits within your card's memory limits.

### 2. GPU-Accelerated Variant Calling (Replaces GATK HaplotypeCaller)
Instead of running GATK HaplotypeCaller on the CPU (`steps/02_variant_calling.sh`), you can call variants at extreme speeds using Clara Parabricks' accelerated haplotype caller:

```bash
pbrun haplotypecaller \
  --ref reference/GCA_000001405.29_GRCh38_no_alt_analysis_set.fna \
  --in-bam out/alignment/SRR622461_merged.bam \
  --recal-file out/alignment/SRR622461_recal.txt \
  --out-vcf out/variants/SRR622461.vcf.gz \
  --num-gpus 1
```

After generating the accelerated BAM and VCF files via Clara Parabricks, you can proceed directly to step `03_qc.sh` and run the rest of the pipeline (QC, annotation, haplogroups, and kinship) as normal.

## Repository Contents & Output Glossary

### Source Code
* `steps/`: Shell scripts corresponding to each stage of the pipeline.
* `lib/`: Helper functions and utilities (includes automated Conda initialization).
* `config.sh`: Main configuration file.
* `run_pipeline.sh`: Master execution script.

### Output Glossary (`out/`)
The `out/` directory contains all generated results, logs, and intermediate files. Below is a guide to what each file and folder means:

* **`haplogroups/`**: Contains the results of the phylogenetic lineage tracing.
  * `haplogroups_summary.tsv`: The main summary table listing the precise paternal (Y-DNA) and maternal (mtDNA) haplogroups for all samples, along with a high-level geographical origin interpretation.
  * `haplogroups.txt`: Raw maternal lineage predictions from HaploGrep3.
  * `classify_Y.nwk`: The reconstructed Newick phylogenetic tree for the Y-chromosome lineages.
* **`kinship/`**: Contains the cohort relatedness and Identity-By-Descent (IBD) results.
  * `kinship_summary.tsv`: The final kinship matrix detailing the exact `PI_HAT` relatedness coefficient between every pair of samples, determining if they are twins, siblings, cousins, or unrelated.
  * `samples_chr9_ibd.genome`: Raw PLINK IBD statistics.
* **`annotation/`**: Clinical variant interpretation files.
  * `chrY_MT_final.hg38_multianno.txt`: The master annotated TSV file containing clinical database cross-references (ClinVar, gnomAD, refGene) for every called variant.
  * `pathogenic.tsv`: Filtered list of variants classified strictly as Pathogenic or Likely Pathogenic.
  * `vus.tsv`: Rare variants classified as Variants of Uncertain Significance.
* **`multiqc/`**: 
  * `multiqc_report.html`: An interactive HTML dashboard aggregating all quality control metrics from FastQC, BWA alignment stats, and QualiMap coverage analysis.
* **`final_bam/` & `bwa/`**: Contains the heavy binary alignment map (`.bam`) files storing the mapped sequencing reads, which are used for visualization in tools like IGV.
* **`chrY_MT_final.vcf.gz`**: The core joint-called Variant Call Format file containing all high-quality SNPs and Indels discovered across the entire cohort.

---

*Built with ❤️ for the bioinformatics community.*

# 🧬 Gut Metagenome QC and Assembly Pipeline

A reproducible **Snakemake** workflow for **metagenomic read quality control, host DNA removal, and assembly** of short-read data (Illumina).  
Designed for **bacterial, archaeal, and fungal** metagenomes from host-associated samples (e.g., mouse or human gut).

---

## 🚀 Features

- **Automatic read download** from SRA/ENA using `fastq-dl`
- **Adapter and quality trimming** using `Trim Galore + FastQC`
- **Host DNA removal** using `BWA + Samtools`
  - Keeps both paired and **singleton** reads to retain low-abundance sequences
- **Compression and cleanup** with `pigz`
- **Assembly** using `metaSPAdes`
  - Optional custom `-k` values
  - Supports unpaired reads via `--s1`
- Ready for large-scale batch runs using SLURM arrays or HPC clusters

---

## 🧩 Workflow Overview

```mermaid
graph TD
A[Download FASTQ] --> B[Trim & QC]
B --> C[Map to Host Genome]
C --> D[Extract Unmapped Reads]
D --> E[Compress Clean Reads]
E --> F[Assemble with metaSPAdes]

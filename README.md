# Nanopore Metagenomic Analysis Pipeline

## 📌 Overview
This repository contains automated SLURM-based shell scripts for processing and analyzing Oxford Nanopore sequencing data. The pipeline is divided into two main workflows:
1. **Quality Control (QC) & Filtering:** Assesses raw data quality, removes low-quality reads, and generates QC reports.
2. **Assembly & Taxonomic Annotation:** Performs *de novo* metagenomic assembly, consensus polishing, and taxonomic identification against the 16S rRNA database.

## 🛠️ Prerequisites and Dependencies
This pipeline is designed for a High-Performance Computing (HPC) environment running the **SLURM** workload manager. 

The following modules and tools are required:
*   **FastQC** (v0.12.1)
*   **Flye** (v2.9.3)
*   **Medaka** (v2.1.1)
*   **BLAST+** (v2.15.0+)
*   **Singularity** (for containerized tools)
    *   `nanofilt.sif` (NanoFilt)
    *   `gfastats.sif` (Gfastats)

## 🚀 Pipeline Workflow

### Phase 1: Quality Control (`quality_raw_data.sh`)
*   **Pre-filtering QC:** Runs FastQC on the raw reads.
*   **Filtering:** Uses NanoFilt (via Singularity) to remove reads with a Q-score below 7.
*   **Post-filtering QC:** Runs FastQC on the cleaned dataset to validate quality improvements.

### Phase 2: Assembly, Polishing & BLAST (`assembly_MAG_and_identification.sh`)
*   **Assembly:** Uses Flye in `--meta` mode for *de novo* metagenomic assembly.
*   **Polishing:** Polishes the assembly using Medaka. Automatically selects the model based on the pore version (R9.4.1 (old) or R10.4.1 (new)).
*   **Metrics:** Extracts assembly statistics using Gfastats (via Singularity).
*   **Taxonomic Annotation:** Downloads the NCBI 16S rRNA database and runs BLASTn on the polished consensus sequences to identify microbial taxa.

## 💻 Usage

### 1. Run Quality Control (`quality_raw_data.sh`)
Submit the QC job by passing the raw FASTQ/FASTQ.gz file as an argument:
```bash
sbatch qc_raw_data.sh /path/to/raw_reads.fastq.gz
```
### 2. Run Assembly and BLAST (`assembly_MAG_and_identification.sh`)

This script handles the *de novo* metagenomic assembly, consensus polishing, and taxonomic identification. It requires two arguments to run successfully: the input file and the sequencing chemistry version (to define the correct Medaka model).

```bash
sbatch assembly_MAG_and_identification.sh <path_to_filtered_fastq.gz> <medaka_version: old|new>
```

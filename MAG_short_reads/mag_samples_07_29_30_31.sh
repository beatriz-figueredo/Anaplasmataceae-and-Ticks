#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=mag_maculatum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=700gb
#SBATCH --output=mag_capys.txt
#SBATCH --error=mag_capys.err
#SBATCH --time=3-00:00:00

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# 1. Load required modules
module load nf-core/3.3.2
module load singularity

# 2. Tell Nextflow where to keep downloaded containers
mkdir -p $HOME/.nextflow/singularity_cache

# Adicionado o 'export' aqui:
export NXF_SINGULARITY_CACHEDIR=$HOME/.nextflow/singularity_cache

# 3. Execução do pipeline
nextflow run nf-core/mag \
    -r 5.4.2 \
    -profile singularity \
    -c custom.config \
    --input samplesheet_mag_capys.csv \
    --outdir results_mag_capys \
    --reads_minlength 50 \
    --fastp_qualified_quality 20 \
    --fastp_trim_polyg \
    --host_fasta /projects/zvbd_lab/Beatriz/MAG_scripts_short_reads/host_genomes/GCA_015741225.1_Hydrochoerus_hydrochaeris_HiC_genomic.fa.gz \
    --bin_domain_classification \
    --prokka_with_compliance \
    --prokka_compliance_centre "ZVBD" \
    --gtdbtk_min_completeness 30 \
    --cat_db /projects/zvbd_lab/Beatriz/MAG_scripts_short_reads/databases/20231120_CAT_gtdb \
    --cat_classify_unbinned

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

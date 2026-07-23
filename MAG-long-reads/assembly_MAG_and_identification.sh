#!/bin/bash

#SBATCH --job-name=bats_assembly_blast
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=%x.job%j.out
#SBATCH --error=%x.job%j.err

# --- CONFIGURATION ---
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Use: sbatch $0 <path_fastq.gz> <medaka_version: old|new>"
    exit 1
fi

INPUT_FILE_FULL_PATH=$(realpath "$1") 
MEDAKA_VERSION_FLAG="$2"
THREADS=20

INPUT_FILENAME=$(basename "$INPUT_FILE_FULL_PATH")
INPUT_DIR=$(dirname "$INPUT_FILE_FULL_PATH")
LOCAL_BASE_DIR=$(dirname "$INPUT_DIR") 

ASSEMBLY_OUTPUT_DIR="$LOCAL_BASE_DIR/assembly_metagenome"
MEDAKA_OUTPUT_DIR="$ASSEMBLY_OUTPUT_DIR/Flye-Medaka"
BLAST_DB_LOCAL_DIR="$LOCAL_BASE_DIR/16S_ribosomal_RNA"
GFA_OUTPUT_PATH="$LOCAL_BASE_DIR/consensus.gfa"
GFASTATS_SIF="$HOME/singularity-containers/gfastats.sif"

mkdir -p "$ASSEMBLY_OUTPUT_DIR" "$MEDAKA_OUTPUT_DIR" "$BLAST_DB_LOCAL_DIR"

# --- Define Medaka Model based on flag ---
if [ "$MEDAKA_VERSION_FLAG" == "old" ]; then
    MEDAKA_MODEL="r941_min_fast_g303"
else
    MEDAKA_MODEL="r1041_e82_400bps_hac_v5.2.0"
fi

echo "================================================="
echo "Starting Pipeline: Assembly -> Polish -> BLAST"
echo "Sample: $INPUT_FILENAME"
echo "================================================="

# --- [1/5] Flye Assembly ---
# Flye often requires anaconda3, which conflicts with miniforge
module purge
module load flye/2.9.3
echo -e "\n[1/5] Running Flye assembly ..."
flye --nano-raw "$INPUT_FILE_FULL_PATH" \
    --out-dir "$ASSEMBLY_OUTPUT_DIR" \
    --threads "$THREADS" \
    --iterations 5 \
    --meta \
    --min-overlap 8000 \
    --genome-size 200000

ASSEMBLY_FASTA="$ASSEMBLY_OUTPUT_DIR/assembly.fasta"

# --- [2/5] Medaka Polishing ---
# Medaka requires miniforge; purging here prevents anaconda3 conflicts
module purge
module load medaka/2.1.1
echo -e "\n[2/5] Running Medaka Polishing..."
medaka_consensus \
    -i "$INPUT_FILE_FULL_PATH" \
    -d "$ASSEMBLY_FASTA" \
    -m "$MEDAKA_MODEL" \
    -o "$MEDAKA_OUTPUT_DIR" \
    -t "$THREADS"

CONSENSUS_FASTA="$MEDAKA_OUTPUT_DIR/consensus.fasta"

# --- [3/5] Gfastats Visualization (via Singularity) ---
module purge
module load singularity/4.1.4
echo -e "\n[3/5] Running gfastats via Singularity..."
singularity run "$GFASTATS_SIF" \
    gfastats "$CONSENSUS_FASTA" -o "$GFA_OUTPUT_PATH"

# --- [4/5 & 5/5] BLAST Setup and Search ---
module load blast/2.15.0+

echo -e "\n[4/5] Preparing 16S rRNA Database..."
if [ ! -f "$BLAST_DB_LOCAL_DIR/16S_ribosomal_RNA.nsq" ]; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz -P "$BLAST_DB_LOCAL_DIR"
    tar -xzf "$BLAST_DB_LOCAL_DIR/16S_ribosomal_RNA.tar.gz" -C "$BLAST_DB_LOCAL_DIR/"
    rm "$BLAST_DB_LOCAL_DIR/16S_ribosomal_RNA.tar.gz"
fi

echo -e "\n[5/5] Executing blastn search..."
blastn \
    -query "$CONSENSUS_FASTA" \
    -db "$BLAST_DB_LOCAL_DIR/16S_ribosomal_RNA" \
    -perc_identity 80 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out "$LOCAL_BASE_DIR/blast16S_assembly.tsv" \
    -outfmt "7 qseqid sseqid stitle pident length slen evalue" \
    -num_threads "$THREADS"

## --- REPORTING RESOURCES ---
echo -e "\n✅ Pipeline completed successfully!"
# Running seff manually after completion is safer if the database has a delay
sh -c "seff $SLURM_JOBID"

#!/bin/bash
#SBATCH --job-name=16S_ONT
#SBATCH --output=%x.job%j.out
#SBATCH --error=%x.job%j.err
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=04:00:00

# --- LOADING MODULES ---
module purge
module load singularity
module load python/3.12.3

# --- CONFIGURATIONS & PATHS ---
SINGULARITY_CONTAINER_HOME=/users/bfiguere/singularity_containers
# Exact names of Singularity images (adjust if extension/name differs)
CHOPPER_IMG="$SINGULARITY_CONTAINER_HOME/chopper.sif"
PORECHOP_IMG="$SINGULARITY_CONTAINER_HOME/porechop.sif"
EMU_IMG="$SINGULARITY_CONTAINER_HOME/emu_3_6_2.sif"

# EMU Database (adjust to your database directory path)
EMU_DB="/projects/zvbd_lab/Beatriz/microbiome_ONT/emu_db"

# --- DECONTAMINATION / FILTER PARAMETERS ---
# Multiplier applied to the negative control's RAW READ COUNTS (not relative
# abundance) before subtracting from each sample's raw counts. E.g. 10 means:
# subtract 10x the read count seen in the negative control from every sample,
# then floor at zero. Increase this for a more conservative (aggressive)
# removal of kitome/contaminant taxa. Raw counts are used instead of relative
# abundance because low-read negative controls can make a single stray read
# look like 100% abundance, which would wrongly wipe out real signal.
SUBTRACTION_FACTOR=10

# Minimum relative abundance to keep after decontamination (0.0010 = 0.10%)
MIN_ABUNDANCE=0.0010

# --- PLOTTING PARAMETERS ---
# Taxonomic level used to group/color the relative abundance stacked bar plot.
GROUP_LEVEL="species"
# How many top taxa (by mean abundance across samples) get their own color;
# the rest are grouped into "Others".
TOP_N_TAXA=15
# Taxonomic order used for the presence/absence plot.
ORDER_OF_INTEREST="Rickettsiales"

# --- INPUT DIRECTORY VALIDATION & ABSOLUTE PATHS ---
if [ -z "$1" ]; then
    echo "Error: Input directory not provided."
    echo "Usage: sbatch $0 /path/to/data_folder [negative_control_barcode] [metadata.tsv]"
    echo "  e.g.: sbatch $0 /path/to/run1 barcode24 /path/to/metadata.tsv"
    exit 1
fi

# 1. Get user input and remove trailing slash (if any)
INPUT_DIR=${1%/}

# 2. Force the system to resolve the real absolute path
INPUT_ABS=$(readlink -f "$INPUT_DIR")

if [ ! -d "$INPUT_ABS" ]; then
    echo "Error: Directory '$INPUT_ABS' does not exist."
    exit 1
fi

# 3. Negative control barcode (optional second argument)
NEG_CONTROL="$2"
if [ -z "$NEG_CONTROL" ]; then
    echo "WARNING: No negative control barcode provided as 2nd argument."
    echo "         Decontamination step will be SKIPPED (only the 0.10% filter will run)."
fi

# 4. Metadata TSV mapping barcode -> sample_name (optional third argument)
# Expected columns (tab-separated): barcode    sample_name
METADATA_TSV="$3"
if [ -n "$METADATA_TSV" ]; then
    METADATA_TSV=$(readlink -f "$METADATA_TSV")
    if [ ! -f "$METADATA_TSV" ]; then
        echo "Error: Metadata file '$METADATA_TSV' does not exist."
        exit 1
    fi
else
    echo "WARNING: No metadata TSV provided as 3rd argument."
    echo "         Plots will use barcode names (e.g. barcode01) instead of sample names."
fi

# --- OUTPUT DIRECTORIES CREATION ---
OUT_BASE="${INPUT_ABS}_results_emu"
CHOPPER_OUT="$OUT_BASE/01_chopper"
PORECHOP_OUT="$OUT_BASE/02_porechop"
EMU_OUT="$OUT_BASE/03_emu"
FILTER_OUT="$OUT_BASE/04_decontam_filtered"
PLOTS_OUT="$OUT_BASE/05_plots"
mkdir -p "$CHOPPER_OUT" "$PORECHOP_OUT" "$EMU_OUT" "$FILTER_OUT" "$PLOTS_OUT"

echo "=========================================================="
echo "Starting QC and Taxonomy pipeline"
echo "Input directory:      $INPUT_ABS"
echo "Output directory:     $OUT_BASE"
echo "Threads:              $SLURM_CPUS_PER_TASK"
echo "Negative control:     ${NEG_CONTROL:-<none>}"
echo "Subtraction factor:   ${SUBTRACTION_FACTOR}x"
echo "Min. abundance kept:  ${MIN_ABUNDANCE} (i.e. $(echo "$MIN_ABUNDANCE * 100" | bc)%)"
echo "=========================================================="

# --- LOOP THROUGH BARCODES ---
for BARCODE_PATH in "$INPUT_ABS"/barcode*; do

    if [ -d "$BARCODE_PATH" ]; then

        BARCODE=$(basename "$BARCODE_PATH")
        echo "--> Processing $BARCODE..."

        # 1. CHOPPER
        echo "    [1/3] Filtering and trimming primers with Chopper..."
        zcat "$BARCODE_PATH"/*.fastq.gz | \
        singularity exec "$CHOPPER_IMG" chopper -q 10 -l 1200 --headcrop 45 --tailcrop 45 > "$CHOPPER_OUT/${BARCODE}_chopped.fastq"

        # 2. PORECHOP
        echo "    [2/3] Removing adapters with Porechop..."
        singularity exec "$PORECHOP_IMG" porechop \
            -i "$CHOPPER_OUT/${BARCODE}_chopped.fastq" \
            -o "$PORECHOP_OUT/${BARCODE}_clean.fastq.gz" \
            --threads "$SLURM_CPUS_PER_TASK" \
            > "$PORECHOP_OUT/${BARCODE}_porechop.log" 2>&1

        rm "$CHOPPER_OUT/${BARCODE}_chopped.fastq"

        # 3. EMU
        echo "    [3/3] Classifying with Emu..."
        mkdir -p "$EMU_OUT/${BARCODE}"

        singularity exec \
            -B "$OUT_BASE":"$OUT_BASE" \
            -B "$EMU_DB":"$EMU_DB" \
            -B "$PWD":"$PWD" \
            "$EMU_IMG" emu abundance \
            "$PORECHOP_OUT/${BARCODE}_clean.fastq.gz" \
            --db "$EMU_DB" \
            --output-dir "$EMU_OUT/${BARCODE}" \
            --threads "$SLURM_CPUS_PER_TASK" \
            --min-align-len 1300 \
            --min-pid 96 \
            --keep-counts \
            --keep-read-assignments

        echo "--> $BARCODE finished."
        echo "----------------------------------------------------------"
    fi
done

echo "All barcodes classified. Starting decontamination and abundance filtering..."

# --- 4. DECONTAMINATION (negative control) + 0.10% ABUNDANCE FILTER ---
# This step merges every barcode's Emu rel-abundance.tsv, optionally subtracts
# a multiple of the negative control's abundance from every other sample
# (kitome/contamination removal), then zeroes out anything below
# MIN_ABUNDANCE and renormalizes.

cat > "$OUT_BASE/decontam_filter.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import glob
import os
import sys
import pandas as pd


NON_TAXA_ROWS = {"unmapped", "mapped_filtered", "mapped_unclassified"}


def load_sample(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", dtype={"tax_id": str})
    df = df[~df["tax_id"].isin(NON_TAXA_ROWS)]
    df = df.set_index("tax_id")
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--emu-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--neg-control", default=None)
    p.add_argument("--subtraction-factor", type=float, default=10.0)
    p.add_argument("--min-abundance", type=float, default=0.001)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    barcode_dirs = sorted(
        d for d in glob.glob(os.path.join(args.emu_dir, "barcode*")) if os.path.isdir(d)
    )
    if not barcode_dirs:
        sys.exit("No barcode directories found inside " + args.emu_dir)

    samples = {}
    for bdir in barcode_dirs:
        barcode = os.path.basename(bdir)
        tsv_files = glob.glob(os.path.join(bdir, "*_rel-abundance.tsv"))
        if not tsv_files:
            print(f"WARNING: no rel-abundance.tsv found for {barcode}, skipping.", file=sys.stderr)
            continue
        samples[barcode] = load_sample(tsv_files[0])

    if not samples:
        sys.exit("No sample tables could be loaded.")

    ref_df = next(iter(samples.values()))
    tax_cols = [c for c in ref_df.columns if c.lower() not in ("abundance", "estimated_counts", "estimated counts")]

    annotation = pd.concat([df[tax_cols] for df in samples.values()])
    annotation = annotation.groupby(annotation.index).first()

    # Raw abundance table (as reported by Emu), kept only for reference.
    abund = pd.DataFrame({bc: df["abundance"] for bc, df in samples.items()}).fillna(0.0)
    abund.to_csv(os.path.join(args.out_dir, "raw_merged_abundance.tsv"), sep="\t")

    # Decontamination is done on RAW READ COUNTS ("estimated counts"), not on
    # relative abundance. Relative abundance is counts / total_reads_of_that_
    # sample: when the negative control has very few total reads (typical for
    # blanks), a single stray read can look like 100% abundance, which would
    # make almost any real sample signal get wrongly zeroed out. Raw counts
    # do not have this distortion.
    counts = pd.DataFrame({bc: df["estimated counts"] for bc, df in samples.items()}).fillna(0.0)
    counts.to_csv(os.path.join(args.out_dir, "raw_merged_counts.tsv"), sep="\t")

    corrected_counts = counts.copy()

    if args.neg_control:
        if args.neg_control not in counts.columns:
            sys.exit(
                f"Negative control barcode '{args.neg_control}' not found among "
                f"samples: {list(counts.columns)}"
            )
        nc_counts = counts[args.neg_control]
        print(
            f"Applying proportional subtraction decontamination (on raw read counts) "
            f"using {args.neg_control} (subtraction factor = {args.subtraction_factor}x)"
        )
        for col in counts.columns:
            if col == args.neg_control:
                continue
            corrected_counts[col] = (counts[col] - args.subtraction_factor * nc_counts).clip(lower=0.0)

        corrected_counts.to_csv(os.path.join(args.out_dir, "decontaminated_counts.tsv"), sep="\t")

        report = pd.DataFrame({"neg_control_raw_count": nc_counts})
        report["max_sample_raw_count"] = counts.drop(columns=[args.neg_control]).max(axis=1)
        report = report[report["neg_control_raw_count"] > 0].sort_values(
            "neg_control_raw_count", ascending=False
        )
        report = report.join(annotation, how="left")
        report.to_csv(os.path.join(args.out_dir, "negative_control_taxa_report.tsv"), sep="\t")
    else:
        print("No negative control provided - skipping decontamination step.")

    # Convert corrected counts back to relative abundance per sample.
    relabund = corrected_counts.copy()
    for col in relabund.columns:
        if args.neg_control and col == args.neg_control:
            continue
        total = relabund[col].sum()
        if total > 0:
            relabund[col] = relabund[col] / total

    relabund.to_csv(os.path.join(args.out_dir, "decontaminated_abundance.tsv"), sep="\t")

    # Apply the 0.10% relative abundance filter and renormalize.
    filtered = relabund.copy()
    filtered[filtered < args.min_abundance] = 0.0
    for col in filtered.columns:
        if args.neg_control and col == args.neg_control:
            continue
        total = filtered[col].sum()
        if total > 0:
            filtered[col] = filtered[col] / total

    final = filtered.join(annotation, how="left")
    final.to_csv(os.path.join(args.out_dir, "final_filtered_abundance.tsv"), sep="\t")

    print("Done. Outputs written to", args.out_dir)


if __name__ == "__main__":
    main()
PYEOF

# --- PYTHON ENVIRONMENT SETUP ---
VENV_DIR="$OUT_BASE/env_16S"

# Checks if the virtual environment already exists; if not, creates a new one
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating Python virtual environment at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

# Activates the virtual environment
source "$VENV_DIR/bin/activate"

# Installs (or ensures they are already installed) the dependencies in the virtual environment
# pip is upgraded to run faster and avoid old version warnings
python3 -m pip install --upgrade pip --quiet
pip install --quiet pandas matplotlib numpy

# Since the environment is now active, python will use its libraries in the steps below:

python3 "$OUT_BASE/decontam_filter.py" \
    --emu-dir "$EMU_OUT" \
    --out-dir "$FILTER_OUT" \
    --neg-control "$NEG_CONTROL" \
    --subtraction-factor "$SUBTRACTION_FACTOR" \
    --min-abundance "$MIN_ABUNDANCE"

echo "Decontamination/filtering outputs written to: $FILTER_OUT"
echo "  - raw_merged_abundance.tsv          -> original relative abundance from Emu, untreated"
echo "  - raw_merged_counts.tsv             -> original raw counts from Emu, untreated"
echo "  - decontaminated_counts.tsv         -> counts after negative control subtraction"
echo "  - decontaminated_abundance.tsv      -> recalculated relative abundance after decontamination"
echo "  - negative_control_taxa_report.tsv  -> taxa observed in the negative control (kitome)"
echo "  - final_filtered_abundance.tsv      -> final table (decontamination + 0.10% filter)"

echo "Generating relative abundance and presence/absence plots..."

# --- 5. PLOTS: RELATIVE ABUNDANCE (stacked bar) + PRESENCE/ABSENCE (order of interest) ---
cat > "$OUT_BASE/plot_taxa.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TAXA_ANNOT_COLS = [
    "species", "genus", "family", "order", "class", "phylum",
    "clade", "superkingdom", "subspecies", "species subgroup", "species group",
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="final_filtered_abundance.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("--group-level", default="genus", choices=["species", "genus", "family"])
    p.add_argument("--order-of-interest", default="Rickettsiales")
    p.add_argument("--metadata", default=None,
                    help="TSV with columns 'barcode' and 'sample_name' to relabel samples in the plots")
    p.add_argument("--neg-control", default=None,
                    help="Barcode ID of the negative control; excluded from all plots if provided")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    barcode_cols = [c for c in df.columns if c.startswith("barcode")]
    if not barcode_cols:
        sys.exit("No barcode* columns found in input table.")

    # Exclude the negative control barcode from every plot: it is only a QC
    # reference used earlier for decontamination and does not represent a
    # real biological sample.
    if args.neg_control:
        if args.neg_control in barcode_cols:
            barcode_cols.remove(args.neg_control)
            print(f"Excluding negative control '{args.neg_control}' from plots.")
        else:
            print(f"WARNING: negative control '{args.neg_control}' not found among "
                  f"barcode columns; nothing excluded.", file=sys.stderr)

    if not barcode_cols:
        sys.exit("No sample columns left to plot after excluding the negative control.")

    # Build display labels for the x-axis (sample names instead of raw
    # barcode IDs), if a metadata table was provided.
    display_labels = list(barcode_cols)
    if args.metadata:
        meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
        if "barcode" not in meta.columns or "sample_name" not in meta.columns:
            sys.exit("Metadata TSV must have 'barcode' and 'sample_name' columns.")
        mapping = dict(zip(meta["barcode"], meta["sample_name"]))
        missing = [bc for bc in barcode_cols if bc not in mapping]
        if missing:
            print(f"WARNING: no sample_name found in metadata for: {', '.join(missing)}. "
                  f"Using barcode ID for these.", file=sys.stderr)
        display_labels = [mapping.get(bc, bc) for bc in barcode_cols]

    # ---------- CHART 1: overall relative abundance (stacked bar) ----------
    group_col = args.group_level
    if group_col not in df.columns:
        sys.exit(f"Column '{group_col}' not found in input table.")

    grouped = df.groupby(group_col)[barcode_cols].sum().fillna(0.0)
    mean_abund = grouped.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abund.head(args.top_n).index.tolist()

    plot_df = grouped.loc[top_taxa].copy()
    others = grouped.drop(index=top_taxa).sum(axis=0)
    if others.sum() > 0:
        plot_df.loc["Others"] = others

    plot_df = plot_df[barcode_cols] * 100  # fraction -> percentage

    fig, ax = plt.subplots(figsize=(max(8, len(barcode_cols) * 0.6), 7))
    x_pos = np.arange(len(barcode_cols))
    bottom = np.zeros(len(barcode_cols))
    cmap = plt.get_cmap("tab20")

    for i, taxon in enumerate(plot_df.index):
        values = plot_df.loc[taxon].values
        ax.bar(x_pos, values, bottom=bottom, label=taxon, color=cmap(i % 20))
        bottom += values

    ax.set_ylabel("Relative abundance (%)")
    ax.set_xlabel("Sample")
    ax.set_title(f"Microbial community composition (top {args.top_n} {group_col})")
    ax.set_ylim(0, 100)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(display_labels, rotation=90)
    legend = ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title=group_col.capitalize())
    for text in legend.get_texts():
        if text.get_text() != "Others":
            text.set_fontstyle("italic")
    plt.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "relative_abundance_stacked_bar.png"), dpi=300)
    plt.close(fig)

    # ---------- CHART 2: presence/absence of the order of interest ----------
    if "order" not in df.columns:
        print("WARNING: 'order' column not found, skipping presence/absence chart.", file=sys.stderr)
        print("Plots written to", args.out_dir)
        return

    sub = df[df["order"].str.lower() == args.order_of_interest.lower()].copy()
    if sub.empty:
        print(f"WARNING: no taxa found for order '{args.order_of_interest}'; presence/absence chart skipped.",
              file=sys.stderr)
        print("Plots written to", args.out_dir)
        return

    label_col = "species" if "species" in sub.columns else group_col
    sub = sub.set_index(label_col)[barcode_cols]
    presence = (sub > 0).astype(int)

    fig2, ax2 = plt.subplots(
        figsize=(max(8, len(barcode_cols) * 0.6), max(3, len(presence) * 0.5))
    )
    cmap2 = matplotlib.colors.ListedColormap(["#f0f0f0", "#c0392b"])
    ax2.imshow(presence.values, aspect="auto", cmap=cmap2, vmin=0, vmax=1)

    ax2.set_xticks(range(len(barcode_cols)))
    ax2.set_xticklabels(display_labels, rotation=90)
    ax2.set_yticks(range(len(presence.index)))
    ax2.set_yticklabels(presence.index, fontstyle="italic")
    ax2.set_xlabel("Sample")
    ax2.set_title(f"Presence/absence of {args.order_of_interest} by sample")

    ax2.set_xticks(np.arange(-0.5, len(barcode_cols), 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(presence.index), 1), minor=True)
    ax2.grid(which="minor", color="white", linewidth=1.5)
    ax2.tick_params(which="minor", bottom=False, left=False)

    plt.tight_layout()
    safe_name = args.order_of_interest.lower().replace(" ", "_")
    fig2.savefig(os.path.join(args.out_dir, f"{safe_name}_presence_absence.png"), dpi=300)
    plt.close(fig2)

    print("Plots written to", args.out_dir)


if __name__ == "__main__":
    main()
PYEOF


PLOT_ARGS=(
    --input "$FILTER_OUT/final_filtered_abundance.tsv"
    --out-dir "$PLOTS_OUT"
    --top-n "$TOP_N_TAXA"
    --group-level "$GROUP_LEVEL"
    --order-of-interest "$ORDER_OF_INTEREST"
)
if [ -n "$METADATA_TSV" ]; then
    PLOT_ARGS+=(--metadata "$METADATA_TSV")
fi
if [ -n "$NEG_CONTROL" ]; then
    PLOT_ARGS+=(--neg-control "$NEG_CONTROL")
fi

python3 "$OUT_BASE/plot_taxa.py" "${PLOT_ARGS[@]}"

echo "Plots written to: $PLOTS_OUT"
echo "  - relative_abundance_stacked_bar.png      -> relative abundance (top ${TOP_N_TAXA} ${GROUP_LEVEL}s + Others), negative control excluded"
echo "  - ${ORDER_OF_INTEREST,,}_presence_absence.png -> presence/absence of ${ORDER_OF_INTEREST} by sample, negative control excluded"
echo "Pipeline completed successfully!"
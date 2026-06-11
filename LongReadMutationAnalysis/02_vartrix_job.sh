#!/bin/bash
# =============================================================================
# Vartrix Variant Genotyping — Per-Barcode Long-Read BAM
# =============================================================================
#
# Description:
#   For each demultiplexed PacBio barcode, extracts cell barcodes from the
#   annotated intersection BAM produced by iso-seq.sh, then runs vartrix to
#   generate a per-cell variant genotype matrix across the MRD mutation panel.
#
#   Steps (per barcode):
#     1. Extract all unique cell barcodes (CB tag) from the annotated BAM
#     2. Index the annotated BAM
#     3. Run vartrix to produce a sparse consensus genotype matrix
#        (rows = variants, cols = cell barcodes; values: 0=no call,
#         1=ref/ref, 2=alt/alt, 3=alt/ref)
#
# Usage:
#   sbatch vartrix_job.sh
#   (or: bash vartrix_job.sh for interactive testing)
#
# Inputs (set in USER CONFIGURATION):
#   - Per-barcode annotated intersection BAMs
#     (<barcode>/<barcode>.annotated.intersect.bam), produced by iso-seq.sh
#   - Reference genome FASTA
#   - Variant VCF (MRD mutation panel, GRCh38)
#   - vartrix binary path
#
# Outputs (written to each barcode subdirectory):
#   - all.cell.barcodes.tsv   : unique CB-tagged cell barcodes in this sample
#   - vartrix_consensus.mtx   : sparse genotype matrix (MatrixMarket format)
#   - vartrix_variants.tsv    : variant list corresponding to matrix rows
#
# Dependencies:
#   samtools, vartrix (v1.1.7+)
#   Conda environment: pb-isoseq
#
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institution.edu


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Demultiplexed barcode IDs — one subdirectory per barcode is expected under
# the working directory, each containing <barcode>.annotated.intersect.bam
barcodes=(
    bc2001--bc2001 bc2009--bc2009 bc2017--bc2017 bc2025--bc2025
    bc2033--bc2033 bc2041--bc2041 bc2049--bc2049 bc2057--bc2057
    bc2065--bc2065 bc2073--bc2073 bc2081--bc2081 bc2089--bc2089
    bc2002--bc2002 bc2010--bc2010 bc2018--bc2018 bc2026--bc2026
    bc2034--bc2034 bc2042--bc2042 bc2050--bc2050 bc2058--bc2058
    bc2066--bc2066 bc2074--bc2074 bc2082--bc2082
)

# Reference genome FASTA (must match the genome used for alignment)
ref_genome="/path/to/ref/fasta/genome.fa"

# Variant VCF — MRD mutation panel (GRCh38 coordinates)
mrd_variants="/path/to/ref/MRD_variant_list.GRCh38.vcf"

# Path to vartrix binary
vartrix="/path/to/tools/vartrix-v1.1.7-x86_64-unknown-linux-gnu/vartrix"

# Working directory — must contain one subdirectory per barcode
work_dir="/path/to/Scripts/iso-seq_pilot"

# Threads for vartrix
THREADS=8


# ── Environment setup ─────────────────────────────────────────────────────────

source activate pb-isoseq
cd "${work_dir}"


# =============================================================================
# Main loop: genotype each barcode with vartrix
# =============================================================================

for barcode in "${barcodes[@]}"; do
    echo "Processing: ${barcode}"
    cd "${barcode}"

    # -- Step 1: Extract unique cell barcodes from CB tag --------------------
    # The CB:Z: tag holds the corrected cell barcode assigned by isoseq correct.
    # Barcodes are extracted, sorted, and deduplicated for use as the vartrix
    # cell barcode allowlist.
    samtools view "${barcode}.annotated.intersect.bam" \
        | awk '{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) {print substr($i,6)}}' \
        | sort | uniq > all.cell.barcodes.tsv

    # -- Step 2: Index the annotated BAM ------------------------------------
    samtools index "${barcode}.annotated.intersect.bam"

    # -- Step 3: Run vartrix ------------------------------------------------
    # Produces a sparse consensus genotype matrix (vartrix_consensus.mtx) and
    # the corresponding variant list (vartrix_variants.tsv).
    # Genotype encoding: 0 = no coverage, 1 = ref/ref, 2 = alt/alt, 3 = alt/ref
    "${vartrix}" \
        -v "${mrd_variants}" \
        -b "${barcode}.annotated.intersect.bam" \
        -f "${ref_genome}" \
        -c all.cell.barcodes.tsv \
        -o vartrix_consensus.mtx \
        --threads "${THREADS}" \
        --out-variants vartrix_variants.tsv

    cd ../
    echo "Done: ${barcode}"
done

echo "All barcodes complete."

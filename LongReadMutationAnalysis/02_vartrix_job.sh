#!/bin/bash
# =============================================================================
# Vartrix Variant Genotyping — Dual Hybrid 
# =============================================================================
#
# Description:
#   For each sample processed by 01_iso-seq.sh, extracts cell barcodes from
#   the mapped BAM and runs vartrix in two modes:
#
#     1. Coverage mode — produces separate ALT and REF read-count matrices,
#        enabling downstream coverage-weighted pseudobulk VAF computation
#     2. Consensus mode — produces a per-cell genotype call matrix
#        (0 = no call, 1 = ref/ref, 2 = alt/alt, 3 = alt/ref)
#
#   Both matrices share the same barcode and variant ordering so they can
#   be loaded and compared directly in the downstream R analysis script.
#
# Usage:
#   sbatch 02_vartrix_job.sh
#
# Inputs:
#   - Mapped BAMs from 01_iso-seq.sh (<results_root>/<sample>/<sample>.mapped.bam)
#   - Reference genome FASTA
#   - Variant VCF (MRD mutation panel, GRCh38)
#
# Outputs (per sample, written to <results_root>/<sample>/vartrix_results/):
#   - all.cell.barcodes.tsv                  : unique CB-tagged cell barcodes
#   - vartrix_alt.mtx / vartrix_ref.mtx      : coverage-mode read count matrices
#   - vartrix_10x_variants.tsv               : variant list (coverage mode)
#   - vartrix_10x_consensus.mtx              : consensus genotype matrix
#   - vartrix_10x_variants_consensus.tsv     : variant list (consensus mode)
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
#SBATCH --job-name=vartrix_pipeline
#SBATCH --output=vartrix_pipeline_%j.out

set -eo pipefail


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Reference genome FASTA (must match the genome used for alignment in 01_iso-seq.sh)
ref_genome="/path/to/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Variant VCF — MRD mutation panel (GRCh38 coordinates)
mrd_variants="/path/to/ref/MRD_allFilteredMuts.GRCh38.vcf"

# Path to vartrix binary
vartrix="/path/to/tools/vartrix-v1.1.7-x86_64-unknown-linux-gnu/vartrix"

# Results root directory (must match out_root from 01_iso-seq.sh)
results_root="/path/to/Scripts/dual_hybrid_batch2"

# Sample names — must match subdirectory names created by 01_iso-seq.sh
SAMPLES=("PAVEPS_IDX" "PAVEPS_MRD" "PAREZH_IDX" "PAXYMC_MRD")

# Number of threads
THREADS=8


# ── Conda setup ───────────────────────────────────────────────────────────────

source /home/user/mambaforge/etc/profile.d/conda.sh
conda deactivate
conda activate pb-isoseq


# =============================================================================
# Per-sample vartrix genotyping
# =============================================================================

for sample in "${SAMPLES[@]}"; do

    bam_file="${results_root}/${sample}/${sample}.mapped.bam"

    echo ""
    echo "========================================"
    echo "  Vartrix: ${sample}"
    echo "  BAM:     ${bam_file}"
    echo "========================================"

    if [[ ! -f "${bam_file}" ]]; then
        echo "WARNING: ${bam_file} not found — skipping ${sample}."
        continue
    fi

    vartrix_dir="${results_root}/${sample}/vartrix_results"
    mkdir -p "${vartrix_dir}"
    cd "${vartrix_dir}"

    # Extract all unique cell barcodes from the CB tag
    echo "[${sample}] Extracting cell barcodes..."
    samtools view "${bam_file}" \
        | awk '{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) {print substr($i,6)}}' \
        | sort -u > all.cell.barcodes.tsv

    # Coverage mode: separate ALT and REF read-count matrices
    # Used for coverage-weighted pseudobulk VAF computation downstream
    echo "[${sample}] Running vartrix (coverage mode)..."
    ${vartrix} \
        -v "${mrd_variants}" \
        -b "${bam_file}" \
        -f "${ref_genome}" \
        -s coverage \
        -c all.cell.barcodes.tsv \
        -o vartrix_alt.mtx \
        --threads "${THREADS}" \
        --out-variants vartrix_10x_variants.tsv \
        --ref-matrix vartrix_ref.mtx

    # Consensus mode: genotype call matrix (0/1/2/3)
    # Used for per-cell genotype assignment and Fisher's exact tests downstream
    echo "[${sample}] Running vartrix (consensus mode)..."
    ${vartrix} \
        -v "${mrd_variants}" \
        -b "${bam_file}" \
        -f "${ref_genome}" \
        -s consensus \
        -c all.cell.barcodes.tsv \
        -o vartrix_10x_consensus.mtx \
        --threads "${THREADS}" \
        --out-variants vartrix_10x_variants_consensus.tsv

    echo "[${sample}] Vartrix done."

done

echo ""
echo "========================================"
echo "  All samples complete."
echo "========================================"

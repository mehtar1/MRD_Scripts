#!/bin/bash
# =============================================================================
# PacBio IsoSeq Hybrid Capture Processing Pipeline
# =============================================================================
#
# Description:
#   Processes a single PacBio HiFi BAM file through the full IsoSeq3 single-
#   cell pipeline, then computes coverage statistics over target mutation loci
#   and bait regions using mosdepth.
#
#   Steps:
#     1. cDNA primer removal and read orientation       (lima)
#     2. UMI and cell barcode clipping                  (isoseq tag)
#     3. Poly-A tail removal and concatemer filtering   (isoseq refine)
#     4. Cell barcode correction against allowlist      (isoseq correct)
#     5. UMI-based deduplication                        (isoseq groupdedup)
#     6. Alignment to reference genome                  (pbmm2)
#     7. Gene and exon annotation of aligned reads      (bedtools intersect)
#     8. Alignment QC metrics                           (samtools)
#     9. Gene and loci coverage                         (bedtools coverage)
#    10. Target and bait coverage with mosdepth         (mosdepth)
#
# Usage:
#   sbatch iso-seq.sh
#   (or: bash iso-seq.sh for interactive testing)
#
# Inputs (set in USER CONFIGURATION):
#   - HiFi BAM file from PacBio sequencer
#   - Reference genome FASTA (GRCh38)
#   - IsoSeq primer FASTA
#   - Cell barcode allowlist (reverse-complemented)
#   - Target variant BED/VCF (mutation loci)
#   - Bait BED (hybridization capture probe regions)
#   - Gene annotation BED files (gene-level and exon-level)
#
# Outputs (written to the working directory):
#   - mapped.sorted.bam              : final aligned, sorted BAM
#   - gene.annotated.intersect.bed   : reads annotated to gene bodies
#   - exon.annotated.intersect.bed   : reads annotated to protein-coding exons
#   - gene_coverage.txt              : per-gene coverage (bedtools)
#   - loci_coverage.txt              : per-bait coverage (bedtools)
#   - target_stats.*                 : per-target-locus mosdepth output
#   - bait_stats.*                   : per-bait-region mosdepth output
#   - MRD_allFilteredMuts.GRCh38.bed : sorted target BED derived from VCF
#   - chrom_order.txt                : contig order extracted from BAM header
#
# Dependencies:
#   lima, isoseq3, samtools, pbmm2, bedtools, mosdepth
#   Conda environment: pb-isoseq
#
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --ntasks=4
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institution.edu


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Input HiFi BAM from the PacBio sequencer
bamfile="/path/to/hifi_reads.bam"

# Reference genome (GRCh38, no-alt analysis set)
fasta="/path/to/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# IsoSeq primer sequences (FASTA)
primers="/path/to/primers.fasta"

# Cell barcode allowlist (reverse-complemented for IsoSeq compatibility)
cell_barcodes="/path/to/gex_737K-arc-v1_REVCOMP.txt"

# Target variant VCF (mutation loci of interest; converted to BED below)
target_vcf="/path/to/MRD_allFilteredMuts.GRCh38.vcf"

# Hybridization capture probe regions BED
bait_bed="/path/to/hybrid_probes.bed"

# Gene annotation BED files (generated from GTF)
gene_bed="/path/to/genes/genes.gtf.gene.bed"
exon_bed="/path/to/genes/exons.proteinCoding.gtf.gene.bed"

# Working directory — all outputs written here
work_dir="/path/to/Scripts/dual_hybrid_pilot"

# Number of threads for mosdepth
THREADS=8


# ── Environment setup ─────────────────────────────────────────────────────────

source /home/user/mambaforge/etc/profile.d/conda.sh
conda deactivate
conda activate pb-isoseq

cd "${work_dir}"


# =============================================================================
# STEP 1–2: cDNA primer removal, read orientation, UMI/barcode clipping
# =============================================================================
# lima removes cDNA primers and orients reads; --per-read retains read-level
# tags required by downstream isoseq commands.
# isoseq tag extracts the UMI (12 bp) and cell barcode (16 bp) using the
# T-12U-16B design string (T=transcript, U=UMI, B=barcode).

lima --per-read --isoseq "${bamfile}" "${primers}" output.bam

isoseq tag output.5p--3p.bam flt.bam --design T-12U-16B


# =============================================================================
# STEP 3–4: Poly-A removal, concatemer filtering, barcode correction
# =============================================================================
# isoseq refine removes poly-A tails and chimeric reads (--require-polya
# enforces that only reads with a poly-A signal are kept as full-length).
# isoseq correct matches raw cell barcodes to the allowlist.

isoseq refine flt.bam "${primers}" fltnc.bam --require-polya

isoseq correct -B "${cell_barcodes}" fltnc.bam corrected.bam


# =============================================================================
# STEP 5: UMI-based deduplication
# =============================================================================
# Reads must be sorted by cell barcode (CB tag) before groupdedup.
# groupdedup collapses reads sharing the same CB+UMI into a single molecule.

samtools sort -t CB corrected.bam -o corrected.sorted.bam
isoseq groupdedup corrected.sorted.bam dedup.bam


# =============================================================================
# STEP 6: Alignment to reference genome
# =============================================================================
# pbmm2 aligns deduplicated molecules to GRCh38 using the ISOSEQ preset,
# which is optimised for long, spliced cDNA reads.
# --split-by-sample writes one BAM per sample tag if multiplexed.

pbmm2 align --preset ISOSEQ --sort --split-by-sample \
    dedup.bam "${fasta}" mapped.bam --log-level INFO


# =============================================================================
# STEP 7: Gene and exon annotation via bedtools intersect
# =============================================================================
# -wao reports all features from the BAM (-a) with overlapping annotation
# from the BED (-b), plus the number of overlapping bases. Unannotated reads
# are retained with a '.' in the annotation columns.

bedtools intersect -wao \
    -a mapped.bam \
    -b "${gene_bed}" \
    -bed > gene.annotated.intersect.bed

bedtools intersect -wao \
    -a mapped.bam \
    -b "${exon_bed}" \
    -bed > exon.annotated.intersect.bed


# =============================================================================
# STEP 8: Alignment QC metrics
# =============================================================================

echo "Unique molecules after deduplication:"
samtools view mapped.bam \
    | awk '{print $1}' | sort | uniq | wc -l

echo "Mean mapping quality (MAPQ):"
samtools view mapped.bam \
    | awk '{sum+=$5; n++} END {if(n>0) print sum/n; else print "No mapped reads."}'

printf '\n'


# =============================================================================
# STEP 9: Sort, index, and compute gene/loci coverage
# =============================================================================
# bedtools coverage requires a coordinate-sorted BAM and matching contig order
# between the BAM and BED files. The contig order is extracted from the BAM
# header and used to re-sort both BED files before coverage computation.

samtools sort mapped.bam -o mapped.sorted.bam
samtools index mapped.sorted.bam

# Extract contig order from BAM header for consistent BED sorting
samtools view -H mapped.sorted.bam \
    | awk '$1=="@SQ"{sub("SN:","",$2); print $2}' > chrom_order.txt

# Sort gene BED to match BAM contig order
bedtools sort -i "${gene_bed}" -faidx chrom_order.txt > gene.sorted.bed
mv gene.sorted.bed "${gene_bed}"

# Sort bait BED to match BAM contig order
bedtools sort -i "${bait_bed}" -faidx chrom_order.txt > bait.sorted.bed
mv bait.sorted.bed "${bait_bed}"

# Convert target VCF to 0-based BED (VCF POS is 1-based; BED is 0-based)
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2-1, $2, $3, $4, $5}' \
    "${target_vcf}" > target.bed
bedtools sort -i target.bed -faidx chrom_order.txt > MRD_allFilteredMuts.GRCh38.bed

# Per-gene and per-loci coverage
bedtools coverage \
    -a "${gene_bed}" -b mapped.sorted.bam -sorted > gene_coverage.txt

bedtools coverage \
    -a "${bait_bed}" -b mapped.sorted.bam > loci_coverage.txt


# =============================================================================
# STEP 10: Mosdepth coverage over target and bait regions
# =============================================================================
# mosdepth reports mean depth per interval. target_stats covers the mutation
# loci; bait_stats covers the hybridization probe footprint.
# Output files: <prefix>.regions.bed.gz and <prefix>.mosdepth.summary.txt

mosdepth hs_stats \
    --by MRD_allFilteredMuts.GRCh38.bed \
    --threads "${THREADS}" \
    mapped.sorted.bam

mosdepth \
    --by <(cut -f1-4 MRD_allFilteredMuts.GRCh38.bed) \
    --threads "${THREADS}" \
    target_stats mapped.sorted.bam

mosdepth \
    --by <(cut -f1-4 "${bait_bed}") \
    --threads "${THREADS}" \
    bait_stats mapped.sorted.bam

echo "All steps complete."

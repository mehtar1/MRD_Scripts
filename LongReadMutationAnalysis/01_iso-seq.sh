#!/bin/bash
# =============================================================================
# PacBio IsoSeq Dual Hybrid Capture Processing Pipeline 
# =============================================================================
#
# Description:
#   Processes demultiplexed PacBio HiFi BAMs through the full IsoSeq3 single-
#   cell pipeline, then annotates aligned reads against gene models, computes
#   coverage over target mutation loci and bait regions, and generates per-read
#   summaries linking each molecule to its cell barcode, UMI, gene, and
#   fragment length.
#
#   The pipeline runs in two loops separated by a one-time BED/VCF sorting
#   step (which requires a mapped BAM for contig order):
#
#   Loop 1 — IsoSeq processing (Steps 1–6):
#     1. cDNA primer removal and read orientation       (lima)
#     2. UMI and cell barcode clipping                  (isoseq tag)
#     3. Poly-A tail removal and concatemer filtering   (isoseq refine)
#     4. Cell barcode correction against allowlist      (isoseq correct)
#     5. UMI-based deduplication                        (isoseq groupdedup)
#     6. Alignment to reference genome                  (pbmm2)
#
#   One-time prep — sort reference BED/VCF to match BAM contig order
#
#   Loop 2 — Annotation, metrics, and coverage (Steps 7–13):
#     7.  Gene and exon annotation                      (bedtools intersect)
#     8.  Alignment QC metrics                          (samtools)
#     9.  Sort + index for coverage                     (samtools)
#     10. Coverage across genes, baits, and targets     (bedtools + mosdepth)
#     11. Per-gene read count summary                   (bedtools + awk)
#     12. Per-read annotation summary (molecule → gene) (samtools + bedtools)
#     13. Reads overlapping mutation sites               (bedtools intersect)
#
# Usage:
#   Edit the USER CONFIGURATION block below, then:
#     sbatch 01_iso-seq.sh
#
# Inputs:
#   - Demultiplexed HiFi BAMs (one per barcode, in bam_dir)
#   - Reference genome FASTA (GRCh38)
#   - IsoSeq primer FASTA
#   - Cell barcode allowlist (reverse-complemented)
#   - Target variant VCF (MRD mutation panel)
#   - Bait BED (hybridization capture probe regions)
#   - Gene, exon, and transcript annotation BED files
#
# Outputs (per sample, written to out_root/<sample>/):
#   - <sample>.mapped.bam                               : aligned BAM
#   - <sample>.gene.annotated.intersect.bed             : reads × gene overlaps
#   - <sample>.exon.annotated.intersect.bed             : reads × exon overlaps
#   - <sample>.gene_coverage.txt                        : per-gene coverage
#   - <sample>.loci_coverage.txt                        : per-bait coverage
#   - <sample>.{hs,target,bait}_stats.*                 : mosdepth output
#   - <sample>.annotated.intersect.geneOverlap.txt      : gene read counts
#   - <sample>.annotated.intersect.summary.tsv           : per-read annotation
#   - <sample>.site_overlap_output.tsv                   : mutation site overlaps
#
# Dependencies:
#   lima, isoseq3, samtools, pbmm2, bedtools, mosdepth
#   Conda environment: pb-isoseq
#
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institution.edu
#SBATCH --job-name=isoseq_pipeline
#SBATCH --output=isoseq_pipeline_%j.out

set -eo pipefail


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Reference genome (GRCh38, no-alt analysis set)
fasta="/path/to/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Input HiFi BAM directory and filename prefix
# Each demultiplexed BAM is at: <bam_dir>/<bam_prefix>.<barcode>--<barcode>.bam
bam_dir="/path/to/Data/PacBio_DualHybridBatch2"
bam_prefix="m64431e_260626_204244.hifi_reads"

# IsoSeq primer sequences and cell barcode allowlist
primers="/path/to/primers.fasta"
cell_barcodes="/path/to/gex_737K-arc-v1_REVCOMP.txt"

# Target variant VCF (MRD mutation panel, GRCh38)
target_vcf="/path/to/ref/MRD_allFilteredMuts.GRCh38.vcf"

# Hybridization capture probe regions BED
bait_bed="/path/to/hybrid_probes.bed"

# Gene annotation BED files
gene_bed="/path/to/ref/genes/genes.gtf.gene.bed"
exon_bed="/path/to/ref/genes/exons.proteinCoding.gtf.gene.bed"
transcript_bed="/path/to/ref/genes/genes.gtf.transcript.bed"

# Output root directory — one subdirectory per sample is created here
out_root="/path/to/Scripts/dual_hybrid_batch2"

# Sample map: PacBio barcode → sample name
# Add or modify entries to match your demultiplexing layout
declare -A SAMPLE_MAP=(
    ["bc2068"]="PAVEPS_IDX"
    ["bc2076"]="PAVEPS_MRD"
    ["bc2084"]="PAREZH_IDX"
    ["bc2092"]="PAXYMC_MRD"
)

# Number of threads for mosdepth
THREADS=8


# ── Conda setup ───────────────────────────────────────────────────────────────

source /home/user/mambaforge/etc/profile.d/conda.sh
conda deactivate
conda activate pb-isoseq


# =============================================================================
# Loop 1: IsoSeq processing (Steps 1–6)
# =============================================================================
# For each barcode: primer removal → tag → refine → correct → dedup → align

for barcode in "${!SAMPLE_MAP[@]}"; do

    sample="${SAMPLE_MAP[$barcode]}"
    input_bam="${bam_dir}/${bam_prefix}.${barcode}--${barcode}.bam"
    sample_dir="${out_root}/${sample}"
    mkdir -p "${sample_dir}"

    echo ""
    echo "========================================"
    echo "  Processing: ${sample}  (${barcode})"
    echo "  Input:      ${input_bam}"
    echo "  Output dir: ${sample_dir}"
    echo "========================================"

    if [[ ! -f "${input_bam}" ]]; then
        echo "WARNING: ${input_bam} not found — skipping ${sample}."
        continue
    fi

    cd "${sample_dir}"

    # Step 1: cDNA primer removal and read orientation
    echo "[${sample}] Step 1: lima (primer removal)"
    lima --per-read --isoseq "${input_bam}" "${primers}" "${sample}.lima.bam"

    # Step 2: Clip UMI (12 bp) and cell barcode (16 bp)
    echo "[${sample}] Step 2: isoseq tag (UMI + cell barcode)"
    isoseq tag "${sample}.lima.5p--3p.bam" "${sample}.tagged.bam" --design T-12U-16B

    # Step 3: Remove poly(A) tails and filter concatemers
    echo "[${sample}] Step 3: isoseq refine (poly-A + concatemer removal)"
    isoseq refine "${sample}.tagged.bam" "${primers}" "${sample}.refined.bam" \
        --require-polya

    # Step 4: Correct cell barcodes against allowlist
    echo "[${sample}] Step 4: isoseq correct (barcode correction)"
    isoseq correct -B "${cell_barcodes}" "${sample}.refined.bam" \
        "${sample}.corrected.bam"

    # Step 5: Sort by CB tag and deduplicate by UMI
    echo "[${sample}] Step 5: groupdedup (UMI deduplication)"
    samtools sort -t CB "${sample}.corrected.bam" -o "${sample}.corrected.sorted.bam"
    isoseq groupdedup "${sample}.corrected.sorted.bam" "${sample}.dedup.bam"

    # Step 6: Align to reference genome (ISOSEQ preset for spliced cDNA reads)
    echo "[${sample}] Step 6: pbmm2 align (genome mapping)"
    pbmm2 align --preset ISOSEQ --sort --split-by-sample \
        "${sample}.dedup.bam" "${fasta}" "${sample}.mapped.bam" --log-level INFO

    cd ..

done


# =============================================================================
# One-time prep: Sort reference BED/VCF to match BAM contig order
# =============================================================================
# Requires a mapped BAM to exist (uses the first sample's BAM for header).
# Creates sorted copies of all reference files in a shared prep/ directory.

first_bam="${out_root}/${SAMPLE_MAP[$(echo ${!SAMPLE_MAP[@]} | awk '{print $1}')]}/"
first_bam=$(ls ${out_root}/*/$(ls ${out_root}/ | head -1)/*.mapped.bam 2>/dev/null | head -1)
# Fallback: use a known sample
first_bam="${out_root}/PAVEPS_IDX/PAVEPS_IDX.mapped.bam"

prep_dir="${out_root}/prep"
mkdir -p "${prep_dir}"

echo "========================================"
echo "  Sorting reference BED files (one-time)"
echo "========================================"

# Extract contig order from BAM header
samtools view -H "${first_bam}" | \
    awk '$1=="@SQ" {sub("SN:","",$2); print $2}' > "${prep_dir}/chrom_order.txt"

# Sort each reference file to match the BAM contig order
# (filter to contigs present in the BAM before sorting)
for src_bed in "${gene_bed}" "${bait_bed}" "${exon_bed}"; do
    out_name=$(basename "${src_bed}" .bed).sorted.bed
    awk 'NR==FNR {c[$1]; next} $1 in c' "${prep_dir}/chrom_order.txt" "${src_bed}" \
        | bedtools sort -i - -faidx "${prep_dir}/chrom_order.txt" \
        > "${prep_dir}/${out_name}"
done

# Convert target VCF to 0-based BED and sort
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2-1, $2, $3, $4, $5}' "${target_vcf}" \
    | awk 'NR==FNR {c[$1]; next} $1 in c' "${prep_dir}/chrom_order.txt" - \
    | bedtools sort -i - -faidx "${prep_dir}/chrom_order.txt" \
    > "${prep_dir}/targets.sorted.bed"

# Simplified gene BED for per-read annotation (chrom, start, end, gene_name)
cut -f1-3,8 "${gene_bed}" > "${prep_dir}/genes.simple.bed"

# Assign sorted file paths for Loop 2
sorted_gene_bed="${prep_dir}/genes.gtf.gene.sorted.bed"
sorted_bait_bed="${prep_dir}/hybrid_probes.sorted.bed"
sorted_exon_bed="${prep_dir}/exons.proteinCoding.gtf.gene.sorted.bed"
sorted_target_bed="${prep_dir}/targets.sorted.bed"
simple_gene_bed="${prep_dir}/genes.simple.bed"


# =============================================================================
# Loop 2: Annotation, metrics, and coverage (Steps 7–13)
# =============================================================================

for barcode in "${!SAMPLE_MAP[@]}"; do

    sample="${SAMPLE_MAP[$barcode]}"
    sample_dir="${out_root}/${sample}"
    cd "${sample_dir}"

    echo ""
    echo "========================================"
    echo "  Post-processing: ${sample}"
    echo "========================================"

    # Step 7: Gene and exon annotation via bedtools intersect
    echo "[${sample}] Step 7: bedtools intersect (gene/exon annotation)"
    bedtools intersect -wao -a "${sample}.mapped.bam" -b "${sorted_gene_bed}" -bed \
        > "${sample}.gene.annotated.intersect.bed"
    bedtools intersect -wao -a "${sample}.mapped.bam" -b "${sorted_exon_bed}" -bed \
        > "${sample}.exon.annotated.intersect.bed"

    # Step 8: Alignment QC metrics
    echo "[${sample}] Step 8: Alignment metrics"
    echo "  Unique molecules after dedup:"
    samtools view "${sample}.mapped.bam" | awk '{print $1}' | sort -u | wc -l
    echo "  Mean MAPQ:"
    samtools view "${sample}.mapped.bam" \
        | awk '{sum+=$5; n++} END {if(n>0) print sum/n; else print "No mapped reads."}'

    # Step 9: Sort + index for coverage tools
    echo "[${sample}] Step 9: Sort and index"
    samtools sort "${sample}.mapped.bam" -o "${sample}.mapped.sorted.bam"
    samtools index "${sample}.mapped.sorted.bam"

    # Step 10: Coverage analysis (bedtools + mosdepth)
    echo "[${sample}] Step 10: Coverage analysis"
    bedtools coverage -a "${sorted_gene_bed}" -b "${sample}.mapped.sorted.bam" \
        > "${sample}.gene_coverage.txt"
    bedtools coverage -a "${sorted_bait_bed}" -b "${sample}.mapped.sorted.bam" \
        > "${sample}.loci_coverage.txt"

    mosdepth --by "${sorted_target_bed}" --threads "${THREADS}" \
        "${sample}.hs_stats" "${sample}.mapped.sorted.bam"
    mosdepth --by <(cut -f1-4 "${sorted_target_bed}") --threads "${THREADS}" \
        "${sample}.target_stats" "${sample}.mapped.sorted.bam"
    mosdepth --by <(cut -f1-4 "${sorted_bait_bed}") --threads "${THREADS}" \
        "${sample}.bait_stats" "${sample}.mapped.sorted.bam"

    # Step 11: Per-gene read count summary
    # Intersects mapped reads with transcript BED, keeps the best gene overlap
    # per read, then counts unique reads per gene.
    echo "[${sample}] Step 11: Gene overlap summary"
    bedtools intersect \
        -a "${sample}.mapped.bam" \
        -b "${transcript_bed}" \
        -bed -wao | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$13,$14,$15,$20,$21,$22}' | \
        sort -k4,4 -k1,1 -k2,2n -k3,3n -k10,10nr | \
        awk '$4!=p{if(p!="") print l; p=$4; l=$0; c=1; next} {c++} c==1{l=$0} END{if(p) print l}' \
        > "${sample}.output.bed"

    cut -f8 "${sample}.output.bed" | sort | uniq -c | sort -nr \
        > "${sample}.annotated.intersect.geneOverlap.txt"

    # Step 12: Per-read annotation summary
    # Extracts molecule ID, cell barcode (CB), UMI (XM), and fragment length
    # from the mapped BAM, then intersects with gene BED to assign each read
    # to its best-overlapping gene. Output: molecule, CB, UMI, gene, overlap, frag_len
    echo "[${sample}] Step 12: Per-read annotation summary"
    samtools view "${sample}.mapped.bam" | \
    awk '
    BEGIN {OFS="\t"}
    NF>=11 && $3!="*" && $1~/molecule/ {
      r=$1;
      c=$3; s=$4-1; e=s+length($10);
      frag_len=e-s;
      cb=""; umi="";
      for(i=11;i<=NF;i++) {
        if($i~/^CB:Z:/) {split($i,t,":"); cb=t[3]}
        if($i~/^XM:Z:/) {split($i,t,":"); umi=t[3]}
      }
      print c,s,e,r,0,"+",cb,umi,frag_len
    }' > "${sample}.reads.bed"

    bedtools intersect \
        -a "${sample}.reads.bed" \
        -b "${simple_gene_bed}" \
        -bed -wao | \
    awk '
    BEGIN {OFS="\t"; print "molecule\tCB\tUMI\tgene\toverlap\tfrag_len"}
    {
      r=$4; cb=$7; umi=$8; g=$13; overlap=$14; frag_len=$9;
      if(overlap==0) next;
      key=r"\t"cb"\t"umi"\t"g;
      if(overlap>(maxo[key]?maxo[key]:0)) {
        maxo[key]=overlap;
        best[key]=g;
        best_frag[key]=frag_len;
      }
    }
    END {
      for(k in maxo) {
        split(k,f,"\t");
        print f[1],f[2],f[3],best[k],maxo[k],best_frag[k]
      }
    }' > "${sample}.annotated.intersect.summary.tsv"

    # Step 13: Identify reads overlapping known mutation sites
    echo "[${sample}] Step 13: Mutation site overlap"
    bedtools intersect \
        -a "${sample}.mapped.bam" \
        -b "${target_vcf}" \
        -wa -wb -bed | \
        sort -k1,1 -k2,2n | \
        uniq -f 5 | \
        cut -f1-20 \
        > "${sample}.site_overlap_output.tsv"

    echo "[${sample}] Done."
    echo ""

done

echo "========================================"
echo "  All samples complete."
echo "========================================"

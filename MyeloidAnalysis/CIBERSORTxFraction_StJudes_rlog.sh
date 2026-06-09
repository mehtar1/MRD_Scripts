# Deconvolution of COG cohort using CIBERSORTx script 

module load R/4.2.3.lua

module load singularity

# Set path to Deconvolution directory (must contain subfolders data/ and outdir/)
datadir=Deconvolution/Projected_Ctype_Decon/StJudes_rlog
# Name of single cell reference file (counts matrix) within data/
ref=Deconvolution/Ctype_Decon/Tables/CtypeAggregate_AllSubtypesBalanced_IDX_RC_NewAnno_SplitEL_MylCombined_1k_111225.txt
# Name of bulk RNA mixture file within data/
mixture=panALL_COGsubset_rlog_111125.txt
# Prefix to save result in Deconvolution directory so that subsequent runs don't overwrite
savePre=StJudes_COG_Fraction0.25_filter_RC_SplitEL_MylComb_ComBatRef

rm -r Deconvolution/Projected_Ctype_Decon/StJudes_rlog/outdir/*

singularity exec --bind /cm/shared --bind $datadir/data:/src/data --bind $datadir/outdir:/src/outdir --bind `pwd` fractions_latest.sif /src/CIBERSORTxFractions --mixture $mixture --rmbatchSmode TRUE --verbose TRUE --perm 500 --refsample $ref --single_cell TRUE --fraction 0.25 --filter FALSE

cp $datadir/outdir/CIBERSORTx_Adjusted.txt $datadir/$savePre\_CIBERSORTx_Adjusted.txt


#!/bin/bash
#SBATCH -J loom
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=n01

source activate mosaic
export LD_LIBRARY_PATH=/share/home/qinxiaohuan/.conda/envs/mosaic/lib:$LD_LIBRARY_PATH


input_dir="/PATH/combined_loom"
output_dir="/PATH/combined_process_loom"
whitelist_dir="/PATH/whitelist_combined"

mkdir "${output_dir}"

for file in ${input_dir}/*.loom; do

    fname=$(basename "$file" .loom)
    fname2=${fname%.combined}

    mkdir "${output_dir}/${fname}"

    python /PATH/preprocess_loom_bc.py \
        -i "${input_dir}/${fname}.loom" \
        -o "${output_dir}/${fname}" \
        --SNP /PATH/1000G_hg19.vcf \
        --ref 37 \
        --region gene \
        --panel /PATH/compass_amplicon_panel.csv \
        --whitelist "${whitelist_dir}/${fname}_whitelist.csv"

done
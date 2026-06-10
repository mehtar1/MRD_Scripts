#!/bin/bash
#SBATCH -J compass
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=n02

source activate mosaic
export LD_LIBRARY_PATH=/share/home/qinxiaohuan/.conda/envs/mosaic/lib:$LD_LIBRARY_PATH

loom_dir=/PATH/combined_loom
indir=/PATH/combined_process_loom
outdir=/PATH/compass_out


patient_info="/PATH/patient_info.tsv"

for file in ${loom_dir}/*.loom; do

    fname=$(basename "$file" .loom)
    fname2=${fname%.combined}

    mkdir -p "$outdir/${fname}"

    patient_sex=$(awk -v target="$fname2" 'BEGIN{FS=OFS="\t"} NR>1 && $2==target {print $10}' "$patient_info")

    if [[ -z "$patient_sex" ]]; then
        echo "Warning: No predicted_sex found for $fname2 in patient_info.tsv"
        continue
    fi

    /PATH/software/COMPASS/COMPASS \
        -i "$indir/${fname}/${fname}" \
        -o "$outdir/${fname}/${fname}" \
        --nchains 4 \
        --chainlength 10000 \
        --CNA 0 \
        --sex "$patient_sex"

    dot -Tpdf -o "$outdir/${fname}/${fname}_tree.pdf" "$outdir/${fname}/${fname}_tree.gv"

done
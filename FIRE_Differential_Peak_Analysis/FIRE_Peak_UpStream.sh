#!/bin/bash
#sBATCH --job-name=test
#SBATCH --partition=gpuh202t112c
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --error=test_%j.err
#SBATCH --output=test_%j.out

# 样本列表，每个样$1本对应两条重复
#samples=(
#"D1L-wt"
#"D24L-wt"
#"D6L-wt"
#"D-wt"
#"L-wt"
#)
source ~/.bashrc
conda init
conda activate chip
samples=($1)

THREADS=5

for ((i=0; i<${#samples[@]}; i++)); do
  for ((j=i+1; j<${#samples[@]}; j++)); do
    s1=${samples[i]}
    s2=${samples[j]}
    pair="${s1}_vs_${s2}"

    echo "============================================"
    echo ">>> Pair: ${pair}"

    beds1=( $(ls ${s1}*peak.bed) )
    beds2=( $(ls ${s2}*peak.bed) )
    allbeds=( "${beds1[@]}" "${beds2[@]}" )

    merged_bed="${pair}_merged_peaks.bed"
    echo "> merge BED: ${allbeds[*]}"

    cat "${allbeds[@]}" \
      | sort -k1,1 -k2,2n \
      | bedtools merge > "${merged_bed}"

    bws1=( $(ls ${s1}*fire.bw ) )
    bws2=( $(ls ${s2}*fire.bw ) )
    allbws=( "${bws1[@]}" "${bws2[@]}" )

    echo "> bigWig files: ${allbws[*]}"

    rawcounts="${pair}_raw_counts.txt"
    npzfile="${pair}.npz"
    echo "> multiBigwigSummary BED-file \\"
    echo "    --bwfiles ${allbws[*]} \\"
    echo "    --BED ${merged_bed} \\"
    echo "    --outRawCounts ${rawcounts} \\"
    echo "    --outFileName ${npzfile}"

    multiBigwigSummary BED-file \
      --bwfiles "${allbws[@]}" \
      --BED "${merged_bed}" \
      --outRawCounts "${rawcounts}" \
      --outFileName "${npzfile}" \
      -p ${THREADS}

    echo ">>> Done pair: ${pair}"
    echo
    conda activate R
    echo "Rscript Peak_Diff.R ${rawcounts} ${s1} ${s2}"
    Rscript Peak_Diff.R ${rawcounts} ${s1} ${s2}
  done
done



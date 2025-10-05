#!/bin/bash
#sBATCH --job-name=test
#SBATCH --partition=tcum256c128Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=52
#SBATCH --error=test_%j.err
#SBATCH --output=test_%j.out

source ~/.bashrc
conda init
conda activate Fiberseq
# 获取日期和样本
date=${1}
sample=${2}
ID=$(basename ${sample} .bam)

# 创建输出目录
output_dir=${date}_${ID}
mkdir -p "${output_dir}"

genome=/home/cxshnl/Data/geneome/B73_V5_Chr/B73_V5_Chr.fa
index=/home/cxshnl/Data/geneome/B73_V5_Chr/B73_V5_Chr.mmi
gene_bed6=/home/cxshnl/Data/geneome/B73_V5_Chr/B73_Chr.bed
len=/home/cxshnl/Data/geneome/B73_V5_Chr/B73_V5_Chr.len
# 执行 pbmm2 align质控--------------------------------------------------------------------------------------------------------
pbmm2 align -j 52 --log-level INFO --sort ${index} ${sample} ${output_dir}/${ID}_mapped.bam 2> ${output_dir}/${ID}.log
samtools view ${sample} | wc -l >> ${output_dir}/${ID}.log

##--------------------------------------------------------------------------------------------------------
total_length=$(samtools stats "${sample}" | grep "total length" | awk '{print $4}')
genome_size=$(seqkit stat ${genome} | awk 'NR==2 {gsub(/,/, "", $5); print $5}')
coverage=$(echo "scale=2; $total_length / $genome_size" | bc)
echo "${sample}: 总碱基数=${total_length} bp" >> ${output_dir}/${ID}.log
echo "测序深度=${coverage}X" >> ${output_dir}/${ID}.log

##----------------------------------------------------------------------------------------------------
samtools view -h -@ 32 ${sample} | awk '{len = length($10); count[len]++} END {for (l in count) print l, count[l]}' > ${output_dir}/${ID}_length.bed
# 执行 m6A 预测

ft predict-m6a -k -t 52 ${output_dir}/${ID}_mapped.bam ${output_dir}/${ID}_m6A_nuc.bam 

ft fire -t 52 -s ${output_dir}/${ID}_m6A_nuc.bam ${output_dir}/${ID}_fire.bam
samtools index ${output_dir}/${ID}_fire.bam

#for i in Chr{1..5}; do
#    samtools view -b ${output_dir}/${ID}_fire.bam ${i} > ${output_dir}/${ID}_m6A_${i}.bam
#    samtools index ${output_dir}/${ID}_m6A_${i}.bam
#    ft fire -t 52 -e --all ${output_dir}/${ID}_m6A_${i}.bam ${output_dir}/${ID}_fire_${i}.bed
#    bedtools sort -i ${output_dir}/${ID}_fire_${i}.bed > ${output_dir}/${ID}_fire_${i}_sorted.bed
    ##awk -v i="$i" '$1 == i' ${output_dir}/${ID}_fire_sorted.bed > ${output_dir}/${ID}_fire_${i}.bed
#    bgzip ${output_dir}/${ID}_fire_${i}_sorted.bed
#    tabix -p bed ${output_dir}/${ID}_fire_${i}_sorted.bed.gz
#done

ft pileup -o ${output_dir}/${ID}_m6A.pileup -m ${output_dir}/${ID}_fire.bam
awk 'NR > 1 { print $1, $2, $3, $9 / $4 }' ${output_dir}/${ID}_m6A.pileup > ${output_dir}/${ID}_m6A.bg
wiggletools mean bin 10 ${output_dir}/${ID}_m6A.bg > ${output_dir}/${ID}_10bp_m6A.bg
bedGraphToBigWig ${output_dir}/${ID}_10bp_m6A.bg ${len}  ${output_dir}/${ID}_10bp_m6A.bw

# 激活 chip 环境
conda activate base
cd ${output_dir}

echo "sample	bam" >> ${ID}.tbl
echo "${ID}	$(realpath "${ID}_m6A_nuc.bam")" >> ${ID}.tbl

echo "ref: ${genome}" >> ${ID}.yaml
echo "ref_name: ${ID}" >> ${ID}.yaml
echo "manifest: $(realpath "${ID}.tbl")" >> ${ID}.yaml
echo "max_t: 52" >> ${ID}.yaml

pixi run fire --configfile ${ID}.yaml

gunzip $(find . -name "*v0.1.1-01-fire-wide-peaks.bed.gz")
High_peak=$(wc -l $(find . -name "*v0.1.1-01-fire-wide-peaks.bed"))
echo "High_Peak : ${High_peak}"  >> ${ID}.log

gunzip $( find . -name "*fire-v0.1.1-wide-peaks.bed.gz")
low_peak=$(wc -l $(find . -name "*fire-v0.1.1-wide-peaks.bed"))
echo "Low_Peak : ${low_peak}"  >> ${ID}.log

cd ./results/${ID}

aligned_bam_to_cpg_scores \
  --bam ${ID}-fire-v0.1.1-filtered.cram \
  --output-prefix ${ID}_CPG \
  --ref ${genome} \
  --threads 52

mv ./trackHub-v0.1.1/bw/all.nucleosome.coverage.bw   ./trackHub-v0.1.1/bw/${ID}_nuc.bw
mv ./trackHub-v0.1.1/bw/all.percent.accessible.bw    ./trackHub-v0.1.1/bw/${ID}_percent.bw
mv ./trackHub-v0.1.1/bw/all.linker.coverage.bw       ./trackHub-v0.1.1/bw/${ID}_linker.bw
mv ./trackHub-v0.1.1/bw/all.fire.coverage.bw         ./trackHub-v0.1.1/bw/${ID}_fire.bw
mv ./trackHub-v0.1.1/bw/hap1.percent.accessible.bw   ./trackHub-v0.1.1/bw/${ID}_hap1.bigwig
mv ./trackHub-v0.1.1/bw/hap2.percent.accessible.bw   ./trackHub-v0.1.1/bw/${ID}_hap2.bigwig
##---------------------------------------------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------------------------------------------
conda activate chip
cd ../..

##----------------------------------------------------------------------------------
bw_files=$(find . -name "*bw*" | sort)
for bw in ${bw_files}; do
  label=$(basename ${bw} | sed 's/.bw//g')
  echo "正在处理: ${label}"

  computeMatrix scale-regions \
    -S ${bw} \
    -R ${gene_bed6} \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    -p 52 \
    -o ${ID}_${label}.gz  # 添加标签到输出文件名

  plotHeatmap -m ${ID}_${label}.gz \
              --plotFileFormat svg \
              --samplesLabel ${label} \
              --heatmapWidth 8 \
              -out ${ID}_${label}.svg

done

echo "所有文件处理完成！"

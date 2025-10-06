input=$1
sample=$2
output_dir="${1}_analysis"


mkdir -p "${output_dir}"


bam_files=$(find ./${input} -name "*${sample}*CLean.bam" | sort)
bam_files_str="${bam_files}"
labels=$(echo ${bam_files} | xargs -n 1 basename | sed 's/_CLean.bam//g' | tr '\n' ' ')



##--------------------------------------------------------------------------------------
bw_files=$(find ./${input} -name "*${sample}*bigwig" | sort)
bw_files_str="${bw_files}"
label1s=$(echo ${bw_files} | xargs -n 1 basename | sed 's/.bigwig//g' | tr '\n' ' ')
echo "BigWig files to process: ${label1s}"


computeMatrix reference-point \
  -S ${bw_files_str} \
  -R Tair_gene_6c.bed \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  -p 35 \
  -o ${output_dir}/${sample}_rp_mat.gz



plotHeatmap -m ${output_dir}/${sample}_rp_mat.gz \
    --samplesLabel ${label1s} \
    --plotFileFormat svg \
    --heatmapWidth 10 \
    -out ${output_dir}/${sample}_rp_heatmap.svg


cd ${output_dir}

macs2_out="${sample}_macs_out"

genome_size="1.2e9" 

bam_files=$(find . -name "*${sample}*.bam")


for bam_file in ${bam_files}; do

    if [[ "${bam_file}" == *"input"* ]]; then
        echo "Skipping ${bam_file} as it contains 'input'"
        continue
    fi


ID1=$(basename ${bam_file} _CLean.bam)
    macs2 callpeak -t ${bam_file} \
                   -c $(find . -name "*${sample}*input*bam") \
                   -g ${genome_size} \
                   -n ${ID1}_Ara \
                   --outdir ${macs2_out} \
                   --format=BAM \
                   --bdg

makeTagDirectory ${ID1}_TagDir ${bam_file}
findPeaks ${ID1}_TagDir -o ${ID1}.peak -gsize 1.2e8 minDist 150 -region

    echo "Peak calling done for ${bam_file}"
done


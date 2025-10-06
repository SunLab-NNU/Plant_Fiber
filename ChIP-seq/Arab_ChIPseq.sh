input="${1}"
hao="${2}"
cd /home/cxshnl/ChIP-Mapping
for sample in $(find . -name "${2}*_R1.fq.gz")
do

 ID=$(basename ${sample} _R1.fq.gz)
 FASTQ1=${sample}
 FASTQ2=${ID}_R2.fq.gz
 THREADS=30
  
 mkdir -p ${input}
 PUT=./${input}
 mkdir -p ${PUT}/${input}_${ID}
 OUTPUT=${PUT}/${input}_${ID}
 mkdir -p ${OUTPUT}/LOG
 LOG=${OUTPUT}/LOG
  
fastp --in1 ${FASTQ1} --in2 ${FASTQ2}  \
--out1 ${OUTPUT}/${ID}_R1.fp.gz \
--out2 ${OUTPUT}/${ID}_R2.fp.gz \
-j ${OUTPUT}/${ID}.json \
-h ${OUTPUT}/${ID}.html \
--thread ${THREADS} \
--length_required 30

Bowtie_Index=/home/cxshnl/Data/geneome/NucTAIR10.fa.genome
CFq_gz1=${OUTPUT}/${ID}_R1.fp.gz
CFq_gz2=${OUTPUT}/${ID}_R2.fp.gz

bowtie2 -q -p ${THREADS} -x ${Bowtie_Index} -1 ${CFq_gz1} -2 ${CFq_gz2} \
    | samtools view -Sb - \
    | samtools sort -@ ${THREADS} -o ${OUTPUT}/${ID}_sorted.bam -

samtools index -@ ${THREADS} ${OUTPUT}/${ID}_sorted.bam

java -jar /home/cxshnl/my_data/picard.jar MarkDuplicates \
I=${OUTPUT}/${ID}_sorted.bam \
O=${OUTPUT}/${ID}_sorted_dedup_reads.bam \
M=${OUTPUT}/${ID}_metrics.txt \
ASSUME_SORTED=TRUE \
CREATE_INDEX=TRUE \

CHROMOSOMES="1 2 3 4 5"

samtools view -@ ${THREADS} -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 ${OUTPUT}/${ID}_sorted_dedup_reads.bam ${CHROMOSOMES} > ${OUTPUT}/${ID}_CLean.bam

samtools index -@ ${THREADS} ${OUTPUT}/${ID}_CLean.bam
samtools idxstats ${OUTPUT}/${ID}_CLean.bam > ${LOG}/${ID}_clean_idxstats.txt &
samtools flagstat ${OUTPUT}/${ID}_CLean.bam > ${LOG}/${ID}_clean_flagstat.txt &

bamCoverage --bam ${OUTPUT}/${ID}_CLean.bam --outFileName ${OUTPUT}/${ID}_bigwig  --outFileFormat bigwig --ignoreDuplicates --binSize 1

done


#!/bin/bash
#sBATCH --job-name=RNA-Seq
#SBATCH --partition=tcum256c128Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --error=DBL_2024.9.16_HS-_%j.err
#SBATCH --output=DBL_2024.9.16_HS-_%j.out

input="${1}"
SAMPLE_NAME="${2}"


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
--length_required 70

CFq_gz1=${OUTPUT}/${ID}_R1.fp.gz
CFq_gz2=${OUTPUT}/${ID}_R2.fp.gz

STAR  --runMode alignReads \
    --genomeDir /home/cxshnl/Data/geneome/Arab.genome.index  \
    --runThreadN ${THREADS}   \
    --readFilesIn ${CFq_gz1} ${CFq_gz2} \
    --outFileNamePrefix ${OUTPUT}/${ID}_STAR \
    --readFilesCommand gunzip -c \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMstrandField intronMotif \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 2  \
    --alignIntronMin 35  \
    --alignIntronMax 2000 \
    --outWigType wiggle \
    --outWigNorm RPM > ${OUTPUT}/${ID}_STAR.log 2>&1
	
samtools index -@ ${THREADS} ${OUTPUT}/${ID}_STARAligned.sortedByCoord.out.bam
BAM="${OUTPUT}/${ID}_STARAligned.sortedByCoord.out.bam"

samtools view -@ ${THREADS} -h ${BAM} |grep -E '^@|NH:i:1' |samtools view -@ ${THREADS} -Sb - > ${OUTPUT}/${ID}_STAR_Uniq.bam
sambamba index -t ${THREADS} ${OUTPUT}/${ID}_STAR_Uniq.bam

samtools index  ${OUTPUT}/${ID}_STAR_Uniq.bam
bamCoverage -b ${OUTPUT}/${ID}_STAR_Uniq.bam \
             -o ${OUTPUT}/${ID}_bigwig.bw \
             --minMappingQuality 30 \
             --smoothLength 50
done

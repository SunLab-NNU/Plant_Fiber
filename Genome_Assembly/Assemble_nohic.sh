#!/bin/bash
#SBATCH --job-name=Assemble_nohic.sh
#SBATCH --partition=tcum512c128Partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=52
#SBATCH --error=test_%j.err
#SBATCH --output=test_%j.out

echo "激活Assembly虚拟环境"
echo "将【要比对的基因组文件】和【Hifi测序的原始bam文件】移动到 ~/Assemble 目录下"
echo "输入 sh this.sh 2XXXXX_species Control.fa sample.bam"

date=$1
Control_FA=$2
query=$3

ID=$(basename ${query} .bam)

# 创建输出目录
output_dir=${date}
mkdir -p "${output_dir}"

samtools fastq -@ 52 ${query} | bgzip > ${ID}.fq.gz
hifiasm -o ./${output_dir}/${ID}.asm -t52 ${ID}.fq.gz
#hifiasm -o ./${output_dir}/${ID}.asm -t 52 -k 63 --lowQ 8- --ctg-n 10 -a 4 -m 10000000 ${ID}.fq.gz
bioconvert gfa2fasta --force -m awk ./${output_dir}/${ID}.asm.bp.p_ctg.gfa ./${output_dir}/${ID}.fasta
python /data/cxshnl/minaconda3/envs/Fiberseq/bin/quast -o ./${output_dir}/${ID}_result ./${output_dir}/${ID}.fasta

cd ${output_dir}

ragtag.py correct ../${Control_FA} ${ID}.fasta -t 52
ragtag.py scaffold ../${Control_FA} ./ragtag_output/ragtag.correct.fasta -t 52
##ragtag.py patch ../${Control_FA} ./ragtag_output/ragtag.scaffold.fasta -t 52

seqkit grep -v -p "Chr0_RagTag" ./ragtag_output/ragtag.scaffold.fasta > filtered_ragtag.fasta
seqkit replace -p "_RagTag" -r "" filtered_ragtag.fasta > renamed_ragtag.fasta
seqkit sort -N renamed_ragtag.fasta > ${ID}_sorted.fa

##ragtag.py patch ./ragtag_output/ragtag.scaffold.fasta scaffolds_FINAL.fasta -t 32
perl /data/cxshnl/Assemble/NGenomeSyn/bin/GetTwoGenomeSyn.pl -InGenomeA ../${Control_FA} -InGenomeB ${ID}_sorted.fa -OutPrefix ${date}

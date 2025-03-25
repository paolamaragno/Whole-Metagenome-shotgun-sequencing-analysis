#!/bin/bash

#PBS -N bowtie2
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/bowtie2_env

outdir="/data/bioinfo/FC/pmaragno/immubac_new/bowtie_out"
workdir="/data/bioinfo/FC/pmaragno/immubac_new/trimmomatic_output"
mkdir ${outdir}

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                bowtie2 -x /data/bioinfo/FC/pmaragno/immubac_new/human_genome/GCF_009914755.1_T2T-CHM13v2.0 -1 ${workdir}/${line}_R1_001_filtered.fastq -2 ${workdir}/${line}_R2_001_filtered.fastq \
                -S ${outdir}/${line}.sam --very-sensitive-local -p 8

        done

conda deactivate
source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/samtools_env

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do
                samtools view -bS ${outdir}/${line}.sam > ${outdir}/${line}.bam
                samtools view -b -f 12 -F 256 ${outdir}/${line}.bam > ${outdir}/${line}.bothunmapped.bam
                samtools sort -n -m 5G -@ 2 ${outdir}/${line}.bothunmapped.bam -o ${outdir}/${line}.bothunmapped.sorted.bam
                samtools fastq ${outdir}/${line}.bothunmapped.sorted.bam -1 >(gzip > ${outdir}/${line}_filtered.final_1.fastq.gz) -2 >(gzip > ${outdir}/${line}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n
                #rm ${s}.sam; rm ${s}.bam; rm ${s}.bothunmapped.bam; rm ${s}.bothunmapped.sorted.bam ### IF YOU WANT TO REMOVE THE INTERMEDIATE FILES
        done

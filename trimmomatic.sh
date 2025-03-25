#!/bin/bash

#PBS -N trimmomatic
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/trimmomatic_env

workdir="/data/bioinfo/FC/pmaragno/immubac_new/fastq"
outdir="/data/bioinfo/FC/pmaragno/immubac_new/trimmomatic_output"

mkdir ${outdir}

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                trimmomatic PE -threads 28 -phred33 -trimlog ${outdir}/${line}_trimmomatic.log ${workdir}/${line}_R1_001.fastq ${workdir}/${line}_R2_001.fastq \
                ${outdir}/${line}_R1_001_filtered.fastq ${outdir}/${line}_R1_001_unpaired.fastq  ${outdir}/${line}_R2_001_filtered.fastq ${outdir}/${line}_R2_001_unpaired.fastq \
                ILLUMINACLIP:/data/bioinfo/FC/trimmomatic_env/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 \
                LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

        done

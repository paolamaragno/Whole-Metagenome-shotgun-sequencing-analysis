#!/bin/bash

#PBS -N cat_fastq
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

workdir="/data/bioinfo/FC/pmaragno/immubac_new/bowtie_out"

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                cat ${workdir}/${line}_filtered.final_1.fastq.gz ${workdir}/${line}_filtered.final_2.fastq.gz > ${workdir}/${line}_filtered.final_R1_R2.fastq.gz
        done

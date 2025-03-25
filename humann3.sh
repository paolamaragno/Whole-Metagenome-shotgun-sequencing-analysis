#!/bin/bash

#PBS -N humann3
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=48:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/metaphlan_env

outdir="/data/bioinfo/FC/pmaragno/immubac_new/humann_out"
workdir="/data/bioinfo/FC/pmaragno/immubac_new/bowtie_out"
metaphlan_out="/data/bioinfo/FC/pmaragno/immubac_new/metaphlan_out"
mkdir ${outdir}

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                humann -i ${workdir}/${line}_filtered.final_R1_R2.fastq.gz --output ${outdir}/${line} --search-mode uniref90 --threads 28 \
                --taxonomic-profile ${metaphlan_out}/${line}_profile.txt --metaphlan /data/bioinfo/FC/metaphlan_env/bin/metaphlan
        done

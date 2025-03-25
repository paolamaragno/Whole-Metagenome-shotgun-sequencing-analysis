#!/bin/bash

#PBS -N metaphlan4.1.1
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/metaphlan_env

db_dir="/data/bioinfo/FC/pmaragno/immubac_new/metaphlan_db"

outdir="/data/bioinfo/FC/pmaragno/immubac_new/metaphlan_out"
workdir="/data/bioinfo/FC/pmaragno/immubac_new/bowtie_out"
mkdir ${outdir}

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                metaphlan ${workdir}/${line}_filtered.final_R1_R2.fastq.gz --input_type fastq  --bowtie2db ${db_dir} --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2out ${outdir}/${line}.bowtie2.bz2 --samout ${outdir}/${line}.sam.bz2 -o ${outdir}/${line}_profile.txt --nproc 28 --read_min_len 36
        done

merge_metaphlan_tables.py ${outdir}/*_profile.txt > ${outdir}/merged_abundance_table.txt

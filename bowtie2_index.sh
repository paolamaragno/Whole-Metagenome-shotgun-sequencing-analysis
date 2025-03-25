#!/bin/bash

#PBS -N bowtie2_index
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/bowtie2_env

outdir="/data/bioinfo/FC/pmaragno/immubac_new/human_genome"

mkdir ${outdir}

bowtie2-build /data/bioinfo/FC/pmaragno/immubac_new/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna ${outdir}/GCF_009914755.1_T2T-CHM13v2.0

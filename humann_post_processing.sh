#!/bin/bash

#PBS -N humann3_post
#PBS -l select=1:ncpus=28:mem=375gb
#PBS -q workq
#PBS -l walltime=48:00:00
#PBS -M paola.maragno98@gmail.com
#PBS -m e

source /data/apps/anaconda3/bin/activate /data/bioinfo/FC/metaphlan_env

outdir="/data/bioinfo/FC/pmaragno/immubac_new/humann_out"

merge_genefamilies="/data/bioinfo/FC/pmaragno/immubac_new/humann_out/merged_genefamilies_KO_renamed"
mkdir ${merge_genefamilies}

merge_pathabundance="/data/bioinfo/FC/pmaragno/immubac_new/humann_out/merged_pathabundance"
mkdir ${merge_pathabundance}

merge_genefamilies_not_renamed="/data/bioinfo/FC/pmaragno/immubac_new/humann_out/merged_genefamilies_KO"
mkdir ${merge_genefamilies_not_renamed}

cat /data/bioinfo/FC/pmaragno/immubac_new/samples.txt | while read line
        do 
                humann_regroup_table -i ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies.tsv -o ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies_KO.tsv --groups uniref90_ko
                humann_rename_table --input ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies_KO.tsv --output ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies_KO_renamed.tsv --names kegg-orthology
                cp ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies_KO_renamed.tsv ${merge_genefamilies}
                cp ${outdir}/${line}/${line}_filtered.final_R1_R2_pathabundance.tsv ${merge_pathabundance}
                cp ${outdir}/${line}/${line}_filtered.final_R1_R2_genefamilies_KO.tsv ${merge_genefamilies_not_renamed}
        done

humann_join_tables --input ${merge_genefamilies} --output ${merge_genefamilies}/all_genefamilies_KO_renamed.tsv --file_name genefamilies_KO_renamed

humann_join_tables --input ${merge_pathabundance} --output ${merge_pathabundance}/all_pathabundance.tsv --file_name pathabundance

humann_join_tables --input ${merge_genefamilies_not_renamed} --output ${merge_genefamilies_not_renamed}/all_genefamilies_KO.tsv --file_name genefamilies_KO

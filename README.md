This repository is a workflow for Whole Metagenome Shotgun sequencing analysis using BioBakery suite.

To perform the analysis run the bash scripts in the following order:
```
qsub trimmomatic.sh
qsub bowtie2_index.sh
qsub bowtie2.sh
qsub cat_pair_ends.sh
qsub metaphlan4.sh
qsub humann3.sh
qsub humann3_post.sh
```

Then run the following command to predict the gut metabolic modules:
```
java -jar /Users/paolamaragno/Downloads/omixer-rpm-1.1.jar -i merged_genefamilies_KO/all_genefamilies_KO_for_omixer.tsv -c 0.66 -d ../GMMs.v1.07.txt
```

Eventually, follow post_processing.R script for all the post-processing analyses in R environment.


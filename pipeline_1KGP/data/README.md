# Example data of 2 genes in 25 samples

The folder provides data snippets to test our pipelines. Counts in the cnt folder can be reproduced from bam files in bam folder using R code ../Docker/get_TReC_ASReC.R. 


- bam - includes bam files for 25 samples from GEUVADIS project in the regions of two genes.

- cnt - the total and allele-specific counts produced from the bam files in folder ```bam```.

- gatkasc - includes output from GATK ASEReadCounter: SNP level allele-specific counts for each sample.


The full genotypes dataset and RNA-seq data were downloaded from arrayexpress database:

- Genotype data: https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/genotypes/

- RNAseq data: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/samples/?query=geuvadis
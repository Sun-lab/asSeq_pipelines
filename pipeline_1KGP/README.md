# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 


This folder contains some data files and codes to run eQTL analysis using the data from 1KGP samples.

- Docker: R code to extract total read count and allele-specific read count

- data: Example data of 2 genes in 25 samples

- datagen: genotype data

- inf: gene information

- runMATRIXEQTL: code to run MATRIXEQTL

- v1_vs_v2: compare asSeq v1 vs. v2 (../asSeq2) where v1 uses standard glm model and v2 uses numerical optimization. 

- test_scr.R: example R code to call ```Docker/get_TReC_ASReC.R``` to generate TReC and ASReC data. 

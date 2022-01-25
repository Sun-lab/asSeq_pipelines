# R code to extract total read count and allele-specific read count 

This R code ```get_TReC_ASReC.R``` was taken from the Docker foler of [gtex_AnVIL repository] (https://github.com/Sun-lab/gtex_AnVIL), with one minor update to change the ```mapqFiler``` from 255 (which was for GTEx bam files) to 20. This is because the RNAseq data in this folder was part of the bam files from Geuvadis Consortium, and the RNA-seq read mapping quality were set in a different fashion. 

```R
ScanBamParam(flag=flag1, what="seq", mapqFilter=20)
```

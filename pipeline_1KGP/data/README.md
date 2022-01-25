

The folder provides data snippets to test our pipelines. Counts in the cnt folder can be reproduced from bam files in bam folder using R code ../Docker/get_TReC_ASReC.R. This code was taken from ```gtex_AnVIL``` repository (https://github.com/Sun-lab/gtex_AnVIL/blob/master/Docker/get_TReC_ASReC.R), with one minor update to change the ```mapqFiler``` from 255 for GTEx bam files to 20. 

```R
ScanBamParam(flag=flag1, what="seq", mapqFilter=20)
```


+ bam - includes bam files for 25 samples from GEUVADIS project in the regions of two genes.

+ cnt - the total and allele-specific counts produced from the bam files in folder ```bam```.

+ gatkasc - includes output from GATK ASEReadCounter: SNP level allele-specific counts for each sample.

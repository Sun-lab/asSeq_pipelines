#  Check the results of dynamic eQTLs associated with CTCF and TP53

## step1_TF.R

Check whether the genes with dynamic eQTLs conditioning on a transcription factor (e.g., CTCF or TP53) are more likley to be one of the targetg genes of the TF. Here we use the target genes defined by Jasper database. The answer is Yes for TP53 and No for CTCF. The insignificant result for CTCF is likely due to the fact that there are only 35 CTCF targe genes. 

## step2_CTCF.R

Check whether genes with dynamic eQTLs conditioning on the expression of CTCF tend to be close to CTCF binding sites. 

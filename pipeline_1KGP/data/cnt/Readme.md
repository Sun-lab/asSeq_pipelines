Information for three example gene processed according to asSeq pipeline, includes
+ TReC_ex.txt - total read counts
+ ASReC_hap1_ex.dat - haplotype 1 allele specific counts
+ ASReC_hap2_ex.dat - haplotype 2 allele specific counts
+ Info_ex.txt - gene information file
+ samples.txt - design matrix information including log(sample depth) and several PCs



Additionally, individual level data including TReC, ASReC_hap1, ASReC_hap2 and ASReC_hapN (conflicting reads) is provided.

Note: In case of running RASQUAL pipeline files Tcnt_ex.dat and samples.dat should be reformatted to be used by RASQUAL: Y.bin, X.bin and K.bin

GE_norm.txt - total counts expression log transformed and quantile normalized for MatrixEQTL analysis and Covariates.txt are samples.txt data transformed to MatrixEQTL format.


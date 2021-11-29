# Detected eQTLs and  Detected dynamic eQTLs using ASE

Each folder contains findings for one particular tissue.

## eQTL results. 

In each tissue-specific folder, there are results for two types of models:

- long model: including all the covariates used in the GTEx study, including: read depth, sex,  genotype PCs, and PEER factors. This is the default model that we used to report the eQTL findings. 

- short model: reduction of the the long model by using only top 5 PEER factors together with all other covariates.

TReCASE package produces the results for joint model, but also for TReC or ASE components of the model separately, thus each model has three output files: the file including the results for the minimum SNP chosen by the full model (using suffix _TReCASE_) as well as for TReC only part of the model (_TReC_) or ASE only part of the model (_ASE_)

Columns include:

gene ids,

MarkerRowID - minimum SNP id

Pos - minimum SNP position

NBod - over-dispersion for Negative-Binomial distribution

BBod - over-dispersion for Negative-Binomial distribution

TReC_b - additive effect size (log fold change) for TReC part of the model

TReC_Pvalue - associated with it p-value

ASE_b - additive effect size (log fold change) for ASE part of the model	

ASE_Pvalue - associated with it p-value

Joint_b	- additive effect size (log fold change) for joint TReCASE model

Joint_Pvalue - associated p-value

trans_Pvalue - test for discordance between TReC and ASE parts of the model

final_Pvalue - final p-value is either the p-value from the joint model or the p-value from the TReC model (if trans_Pvalue is significant)

permp - permutation p-value (for final_Pvalue in TReCASE related file, for TReC_Pvalue for TReC related file and for ASE_Pvalue for ASE related file)

## Dynamic eQTLs are stored in BetaBin folder:

For three contexts of dynamic eQTL or factors of interest (age,  CTCF expression and TP53 expression) we fitted the model searching for changes of eQTL size.

This model is fit only for the eGenes with permutation p-value from the TReCASE model smaller than 0.01 and we only consider the SNP with minimum p-value for each eGene. 

We fit 3 models: 

- nocov - only an association between the factor of interest and allele-specific expression,

- gPC - association between the factor of interest and allele-specific expression including first two gPCs as covariates,

- gPC_PF - association between the factor of interest and allele-specific expression including first two gPCs and first three Peer Factors as covariates,

Columns include:

nm - gene id,

od - over-dispersion (NA if binomial model is used due to low over-dispersion)

int, b.gPC1, b.gPC2, b.PF1, b.PF2, b.PF3, b.PF4, b.PF5 - covariates included in the model

b.cnd - condition of interest (age, CTCF or TP53 in these examples)

e.int, e.od, e.gPC1, e.gPC2, e.PF1, e.PF2, e.PF3, e.PF4, e.PF5, e.cnd - respective errors

pval - likelihood ratio test for the condition of interest,

status - 0 if beta-binomial model was fitted and 1 if binomial model was fitted

qval - qvalue produced from pval

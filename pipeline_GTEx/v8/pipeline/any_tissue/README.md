Prerequisite files in the base folder:
(say, /pine/scr/z/h/zhabotyn/R01/GTEx/v8)

There should be following folders:
1. Reference - contains gencode.v26.GRCh38.genes.gtf, as well as information for CTCF and TP53 sites
2. TReC_ASReC - total and allele-specific expression (in folders named accordingly to the tissue)
3. Annotations/GTEx_Analysis_v8_eQTL_covariates - covariates for each tissue
4. WGS_VCF - directory with phased vcf file

Also, note, that pipeline relies on 2 external packages that need to be preinstalled and their path should be added to the specifications file:
eigenMT and parser4grNAunf - first one to do a fast initial eigenMT estimate of number of effective tests and second to deal with memory intensive processing of vcf files

Once we have the data structure ready, we can run the following pipeline:

#0a. run (1 process):
#at this step we fix suspicious allele-specific expression and fix influential total expression entries for the long model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long' step0_submit_trim.R step0_submit_trim_long_AVO.Rout

#0b. When previous step is finished, run (1 process):
#at this step we fix influential total expression entries for the short model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt short' step0_submit_trim.R step0_submit_trim_long_AVO.Rout

#1. When previous step is finished, run (22 processes)
#at this step we produce preformatted data to be run at a gene level
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step1_submit_preprocSNP.R step1_submit_preprocSNP.Rout

#2a. When previous step is finished, run (#genes processes)
#fitting TReCASE model for each gene, long model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long 5e5' step2_submit_trecaseA.R step2_submit_trecaseA_long.Rout

#2b. As long as there is place in the scheduler run (%genes processes) 
#since schedulers often have something like 30-50k hard limit may need to wait until 2a is done
#fitting TReCASE model for each gene, short model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt short 5e5' step2_submit_trecaseA.R step2_submit_trecaseA_short.Rout

#3. 1 process
#basic comparison of long and short models
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt 5e5' step3_submit_brief_summ.R step3_submit_brief_summ.Rout

#4a. starts 22 processes, that parallelize into #gene processes to
#produce the data to estimate permutation p-value for the long model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long 5e5' step4_submit_MatrixEQTL.R step4_submit_MatrixEQTL_long.Rout

#4b. starts 22 processes, that parallelize into #gene processes to
#same for the short model
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt short 5e5' step4_submit_MatrixEQTL.R step4_submit_MatrixEQTL_short.Rout

#5. collect the results from long and short models (1 process)
#note, here we also reformat allele-specific counts for future analysis of dynamic ASE analysis
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt 5e5' step5_submit_collect.R step5_submit_collect.Rout

#6. #12 marginal glm fitted - gPC1,gPC2, and first ten PEER factors
#for 3 conditions: AGE, CTCF and TP53: a model fit along with gPC1,gPC2 and
#along with gPC1,gPC2, PEER1,...PEER5
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step6_submit_glm.R step6_submit_glm.Rout

#7. summarizing the results in a few plots
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step7_submit_plot.R step7_submit_plot.Rout

#8.
#checking for enrichment for the results step 6
R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step6_submit_glm.R step6_submit_glm.Rout




library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(asSeq)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("four arguments are expected\n", call.=FALSE)
} else {
  bam_file       = args[1]
  gene_anno_file = args[2]
  sam_name       = args[3]
  het_snp_file   = args[4]
}

# setwd("~/research/GitHub/gtex_AnVIL_data/Docker_testing")
# tag = "GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files_"
# sam_name = "GTEX-1117F-0426-SM-5EGHI"
# 
# bam_file = paste0(tag, sam_name, ".Aligned.sortedByCoord.out.patched.md.bam")
# 
# gene_anno_file = "exon_by_genes_gencode.v26.GRCh38.rds"
# gene_anno_file = paste0("~/research/data/_GTEx/v8/Reference/", gene_anno_file)
# het_snp_file   = "GTEX-1117F.txt"

print(sprintf("sam_name: %s", sam_name))

bam_filtered   = paste0(sam_name, "_filtered.bam")
bam_fS         = paste0(sam_name, "_filtered_sorted_byQname")

# ------------------------------------------------------------------------
# counting
# ------------------------------------------------------------------------

ct1 = countBam(bam_file)
print("done with first counting!\n")

# ------------------------------------------------------------------------
# getUnique and filtering
# ------------------------------------------------------------------------

flag1  = scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE,
                     isDuplicate=FALSE, isNotPassingQualityControls=FALSE,
                     isSupplementaryAlignment=FALSE, isProperPair=TRUE)

param1 = ScanBamParam(flag=flag1, what="seq", mapqFilter=20)

filterBam(bam_file, destination=bam_filtered, param=param1)

print("done with filtering!")

# ------------------------------------------------------------------------
# counting again
# ------------------------------------------------------------------------

ct2 = countBam(bam_filtered)
print("done with second counting!\n")

print("the total number of reads/nucleotides before/after filtering:")
print(ct1)
print(ct2)

# ------------------------------------------------------------------------
# calculate total read count (TReC) per gene
# ------------------------------------------------------------------------

genes   = readRDS(gene_anno_file)
bamfile = BamFileList(bam_filtered, yieldSize=1000000)

se = summarizeOverlaps(features=genes, reads=bamfile, mode="Union",
             singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

ct = as.data.frame(assay(se))

print("done with TReC!")

# ------------------------------------------------------------------------
# sort bam files, this is needed for extracting allele-specific reads
# maximum memory is set to be 16G here
# ------------------------------------------------------------------------

date()
sortBam(bam_filtered, bam_fS, byQname=TRUE, maxMemory=3500)
date()

print("done with sortBam!")

# ------------------------------------------------------------------------
# extract allele-specific reads
# ------------------------------------------------------------------------

date()
extractAsReads(input=paste0(bam_fS, ".bam"), snpList=het_snp_file,
               min.avgQ=20,  min.snpQ=20)
date()

print("done with extractAsReads!")

# ------------------------------------------------------------------------
# count allele-specific reads
# ------------------------------------------------------------------------

se1 = summarizeOverlaps(features=genes, reads=paste0(bam_fS, "_hap1.bam"),
                        mode="Union", singleEnd=FALSE, ignore.strand=TRUE,
                        fragments=TRUE)

se2 = summarizeOverlaps(features=genes, reads=paste0(bam_fS, "_hap2.bam"),
                        mode="Union", singleEnd=FALSE, ignore.strand=TRUE,
                        fragments=TRUE)

seN = summarizeOverlaps(features=genes, reads=paste0(bam_fS, "_hapN.bam"),
                        mode="Union", singleEnd=FALSE, ignore.strand=TRUE,
                        fragments=TRUE)

print("done with ASReC!")

ct1 = as.data.frame(assay(se1))
ct2 = as.data.frame(assay(se2))
ctN = as.data.frame(assay(seN))

if(! all(rownames(ct) == rownames(ct1))){
  stop("rownames of ct and ct1 do not match\n")
}

if(! all(rownames(ct1) == rownames(ct2))){
  stop("rownames of ct1 and ct2 do not match\n")
}

if(! all(rownames(ct1) == rownames(ctN))){
  stop("rownames of ct1 and ctN do not match\n")
}

cts = cbind(ct, ct1, ct2, ctN)
dim(cts)
cts[1:2,]

write.table(cts, file = sprintf("%s.trecase.txt", sam_name), 
            quote = FALSE, sep = "\t", eol = "\n")

print("done!")

sessionInfo()
q(save="no")

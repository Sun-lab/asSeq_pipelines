
library(data.table)
library(GenomicRanges)

# --------------------------------------------------------------------
# read in results
# --------------------------------------------------------------------

tissues = list.dirs("../../Results", full.names=FALSE, 
                    recursive=FALSE)

ctcf_all = NULL
genes_all = NULL

for(t1 in tissues){
  ctcf_file = sprintf("%s_ctcf_long_gPC_PF_01.csv", t1)
  ctcf = fread(sprintf("../../Results/%s/BetaBin/%s", t1, ctcf_file))
  n1 = nrow(ctcf)
  genes_all = c(genes_all, ctcf$nm)

  ctcf = ctcf[ctcf$qval < 0.1,.(nm, pval, status, qval)]
  n2 = nrow(ctcf)
  ctcf$tissue = rep(t1, nrow(ctcf))
  
  ctcf_all = rbind(ctcf_all, ctcf)
  cat(sprintf("%s, %d, %d\n", t1, n1, n2))
}

ctcf_all[1:2,]

table(ctcf_all$tissue)

t2 = table(ctcf_all$nm)
table(t2)
t2[t2==2]

t3 = table(genes_all)
table(t3)

genes_all = unique(genes_all)
length(genes_all)

# --------------------------------------------------------------------
# read in gene annotation information
# --------------------------------------------------------------------

gene_file = "gencode.v26.GRCh38.genes_gene_level_anno.txt"
gene_file = file.path("../../Reference", gene_file)
genes = fread(gene_file, data.table=FALSE)
dim(genes)
genes[1:2,]

table(ctcf_all$nm %in% genes$geneId)
table(genes_all %in% genes$geneId)

mat1 = match(unique(ctcf_all$nm), genes$geneId)

genes$ctcf_eQTL = rep(FALSE, nrow(genes))
genes$ctcf_eQTL[mat1] = TRUE
table(genes$ctcf_eQTL)

genes$promoter_start = rep(NA, nrow(genes))
genes$promoter_end   = rep(NA, nrow(genes))

wn = which(genes$strand == "-")
wp = which(genes$strand == "+")

genes$promoter_start[wn] = genes$end[wn]
genes$promoter_end[wn]   = genes$end[wn] + 199

genes$promoter_start[wp] = genes$start[wp] - 199
genes$promoter_end[wp]   = genes$start[wp]

dim(genes)
genes[1:5,]

gr1 = makeGRangesFromDataFrame(genes, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="promoter_start", 
                               end.field="promoter_end")

genes_eQTL = genes[genes$geneId %in% genes_all,]
dim(genes_eQTL)
genes_eQTL[1:5,]

gr0 = makeGRangesFromDataFrame(genes_eQTL, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="promoter_start", 
                               end.field="promoter_end")

# --------------------------------------------------------------------
# read in CTCF binding site information
# --------------------------------------------------------------------

ff1  = "CTCFBSDB_all_exp_sites_Sept12_2012_hg38_loci.bed.gz"

ctcf.bs = fread(file.path("CTCF",ff1))
dim(ctcf.bs)
ctcf.bs[1:5,]

names(ctcf.bs) = c("chr", "start", "end")
table(ctcf.bs$chr)

lens = ctcf.bs$end - ctcf.bs$start + 1
summary(lens)

pdf("CTCF/CTCFBS_len_hist.pdf", width=6, height=4)
par(mar=c(5,4,1,1), bty="n")
hist(log10(lens), xlab="log10(CTCF BS length)", main="", breaks=100)
dev.off()

table(lens < 400)/length(lens)
table(lens < 300)/length(lens)
table(lens < 200)/length(lens)

ctcf.bs$start100 = round(0.5*(ctcf.bs$start + ctcf.bs$end)) - 100
ctcf.bs$end100   = round(0.5*(ctcf.bs$start + ctcf.bs$end)) + 100

gr2 = makeGRangesFromDataFrame(ctcf.bs, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", 
                               end.field="end")

gr3 = makeGRangesFromDataFrame(ctcf.bs, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start100", 
                               end.field="end100")

fun1 <- function(x) sum(width(reduce(x, ignore.strand=T)))
fun1(gr1)

width2 = fun1(gr2)
width3 = fun1(gr3)
width2
width3

prp2 = width2/(3234.83*10^6)
prp3 = width3/(3234.83*10^6)
prp2
prp3

# --------------------------------------------------------------------
# check overlap for all annotated genes
# --------------------------------------------------------------------

mtch2 = findOverlaps(gr1, gr2, select="first")
table(!is.na(mtch2))

mtch3 = findOverlaps(gr1, gr3, select="first")
table(!is.na(mtch3))

table(genes$ctcf_eQTL, !is.na(mtch2))
table(genes$ctcf_eQTL, !is.na(mtch3))

c1 = chisq.test(genes$ctcf_eQTL, !is.na(mtch2))
c2 = chisq.test(genes$ctcf_eQTL, !is.na(mtch3))

c1
c1$expected
c1$observed

c2
c2$expected
c2$observed

# --------------------------------------------------------------------
# check overlap for all genes used in eQTL
# --------------------------------------------------------------------

mtch2 = findOverlaps(gr0, gr2, select="first")
table(!is.na(mtch2))

mtch3 = findOverlaps(gr0, gr3, select="first")
table(!is.na(mtch3))

c1 = chisq.test(genes_eQTL$ctcf_eQTL, !is.na(mtch2))
c2 = chisq.test(genes_eQTL$ctcf_eQTL, !is.na(mtch3))

c1
c1$expected
c1$observed

c2
c2$expected
c2$observed

# --------------------------------------------------------------------
# Compare q-value between genes with or without CTCF binding sites
# --------------------------------------------------------------------

dim(ctcf_all)
ctcf_all[1:2,]

genes_with_CTCF = genes_eQTL$geneId[!is.na(mtch2)]
length(genes_with_CTCF)

ctcf_all$CTCF = rep(0, nrow(ctcf_all))
ctcf_all$CTCF[which(ctcf_all$nm %in% genes_with_CTCF)] = 1
table(ctcf_all$CTCF)


genes_with_CTCF = genes_eQTL$geneId[!is.na(mtch3)]
length(genes_with_CTCF)

ctcf_all$CTCF_200bp = rep(0, nrow(ctcf_all))
ctcf_all$CTCF_200bp[which(ctcf_all$nm %in% genes_with_CTCF)] = 1
table(ctcf_all$CTCF_200bp)

wilcox.test(ctcf_all$qval ~ ctcf_all$CTCF)
wilcox.test(ctcf_all$qval ~ ctcf_all$CTCF_200bp)

# --------------------------------------------------------------------
# Check ASReC data
# --------------------------------------------------------------------

table(ctcf_all$tissue == "Whole_Blood")
genes2check = ctcf$nm[ctcf$tissue == "Whole_Blood"]

fnm = "../ASE_Whole_Blood_counts/Whole_Blood_preprASE_long.csv.gz"
asrec = fread(fnm, header=TRUE)

table(genes2check %in% asrec$V1)

dim(asrec)
asrec[1:2,1:5]

gc()
sessionInfo()
q(save = "no")




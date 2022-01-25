
library(data.table)
library(GenomicRanges)

TF = "TP53"

# --------------------------------------------------------------------
# read in results
# --------------------------------------------------------------------

tissues = list.dirs("../../Results", full.names=FALSE, 
                    recursive=FALSE)

ctcf_all = NULL

for(t1 in tissues){
  ctcf_file = sprintf("%s_ctcf_long_gPC_PF_01.csv", t1)
  ctcf = fread(sprintf("../../Results/%s/BetaBin/%s", t1, ctcf_file))
  n1 = nrow(ctcf)
  
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

table(ctcf_all$status)
tapply(ctcf_all$qval, ctcf_all$status, summary)

# --------------------------------------------------------------------
# read in gene annotation information
# --------------------------------------------------------------------

gene_file = "gencode.v26.GRCh38.genes_gene_level_anno.txt"
gene_file = file.path("../../Reference", gene_file)
genes = fread(gene_file)
dim(genes)
genes[1:2,]

table(ctcf_all$nm %in% genes$geneId)

mat1 = match(unique(ctcf_all$nm), genes$geneId)
ctcf_genes = genes$hgnc_symbol[mat1]
length(ctcf_genes)
table(ctcf_genes == "")

ctcf_genes = ctcf_genes[which(ctcf_genes != "")]
cat(ctcf_genes, sep="\n")

# --------------------------------------------------------------------
# read in CTCF target gene information
# --------------------------------------------------------------------

dir1 = "~/research/data/human/harmonizome/CHEA"
ff1  = "gene_attribute_matrix.txt.gz"
ff1  = file.path(dir1, ff1)

tf.gene = fread(ff1, data.table=FALSE, na.strings = "na")
dim(tf.gene)
tf.gene[1:5,1:5]

table(tf.gene$CTCF)

ctcf_targets = tf.gene[which(tf.gene$CTCF == 1),1]
length(ctcf_targets)

table(ctcf_genes %in% ctcf_targets)

names(ctcf.bs) = c("chr", "start", "end")
table(ctcf.bs$chr)

lens = ctcf.bs$end - ctcf.bs$start + 1
summary(lens)

pdf("../../Reference/CTCF/CTCFBS_len_hist.pdf", width=6, height=4)
par(mar=c(5,4,1,1), bty="n")
hist(log10(lens), xlab="log10(CTCF BS length)", main="", breaks=100)
abline(v=log10(200))
dev.off()

table(lens < 400)/length(lens)
table(lens < 300)/length(lens)
table(lens < 200)/length(lens)

gr2 = makeGRangesFromDataFrame(ctcf.bs, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", 
                               end.field="end")

gr3 = makeGRangesFromDataFrame(ctcf.bs[which(lens < 200),], 
                               ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", 
                               end.field="end")

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

mtch2 = findOverlaps(gr1, gr2, select="first")
table(!is.na(mtch2))

mtch3 = findOverlaps(gr1, gr3, select="first")
table(!is.na(mtch3))

table(ctcf$qval < 0.1, !is.na(mtch2))
table(ctcf$qval < 0.1, !is.na(mtch3))

table(ctcf$qval < 0.25, !is.na(mtch2))
table(ctcf$qval < 0.25, !is.na(mtch3))

chisq.test(ctcf$qval < 0.1, !is.na(mtch2))
chisq.test(ctcf$qval < 0.1, !is.na(mtch3))

gc()
sessionInfo()
q(save = "no")




library(data.table)

TF = "CTCF"
tf = tolower(TF)

# --------------------------------------------------------------------
# read in results
# --------------------------------------------------------------------

tissues = list.dirs("../../Results", full.names=FALSE, 
                    recursive=FALSE)

tf_all = NULL

for(t1 in tissues){
  tf_file = sprintf("%s_%s_long_gPC_PF_01.csv", t1, tf)
  dyn_eqtl = fread(sprintf("../../Results/%s/BetaBin/%s", t1, tf_file))
  n1 = nrow(dyn_eqtl)
  
  dyn_eqtl = dyn_eqtl[dyn_eqtl$qval < 0.1,.(nm, pval, status, qval)]
  n2 = nrow(dyn_eqtl)
  dyn_eqtl$tissue = rep(t1, nrow(dyn_eqtl))
  
  tf_all = rbind(tf_all, dyn_eqtl)
  cat(sprintf("%s, %d, %d\n", t1, n1, n2))
}

tf_all[1:2,]

table(tf_all$tissue)

t2 = table(tf_all$nm)
table(t2)
t2[t2==2]

# --------------------------------------------------------------------
# read in gene annotation information
# --------------------------------------------------------------------

gene_file = "gencode.v26.GRCh38.genes_gene_level_anno.txt"
gene_file = file.path("../../Reference", gene_file)
genes = fread(gene_file)
dim(genes)
genes[1:2,]

table(tf_all$nm %in% genes$geneId)

mat1 = match(unique(tf_all$nm), genes$geneId)
tf_genes = genes$hgnc_symbol[mat1]
length(tf_genes)
table(tf_genes == "")

tf_genes = tf_genes[which(tf_genes != "")]
cat(tf_genes, sep="\n")

# --------------------------------------------------------------------
# read in CTCF target gene information
# --------------------------------------------------------------------

dir1 = "~/research/data/human/harmonizome/Jaspar\ PWMs"
ff1  = "gene_attribute_matrix.txt.gz"
ff1  = file.path(dir1, ff1)

target_anno = fread(ff1, data.table=FALSE, na.strings = "na")
dim(target_anno)
target_anno[1:5,1:5]

table(target_anno$TP53)
table(target_anno$CTCF)

tf_targets = target_anno[which(target_anno[[TF]] == 1),1]
length(tf_targets)

all_targets = target_anno[-(1:3),1]

table(tf_genes %in% tf_targets)
table(all_targets %in% tf_genes, all_targets %in% tf_targets)
chisq.test(all_targets %in% tf_genes, all_targets %in% tf_targets)

gc()
sessionInfo()
q(save = "no")



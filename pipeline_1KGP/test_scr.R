library("GenomicFeatures")

#
#create gene annotation
gtfFile = "inf/gencode.v26.GRCh38.genes.gtf"

genes = read.table(gtfFile, sep="\t", as.is=T)
genes = genes[genes[,3]=="gene",]
genes[1:2,]

gnkp = c("ENSG00000008128.22", "ENSG00000008130.15", "ENSG00000189339.11")
get_block = function(x, split=";", block=1){
  unlist(strsplit(x, split=split))[block]
}
geneids = sapply(sapply(genes[,9], get_block), get_block, split=" ", block=2)
names(geneids) = NULL
genes = genes[geneids %in% gnkp,]
geneids = geneids[geneids %in% gnkp]
gene_info = data.frame(id=geneids, chr=genes[,1], start=genes[,4], end=genes[,5])
write.table(gene_info, "cnt/Info_ex.txt", quote=F, sep="\t", row.names=F)

path = "https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf/"

txdb = makeTxDbFromGFF(file=gtfFile, format="gtf",
  dataSource=paste(path, gtfFile, sep=""),
  organism="Homo sapiens")

saveDb(txdb, file="gencode.v26.GRCh38.genes.sqlite")

seqlevels(txdb)
columns(txdb)
keytypes(txdb)

genes = exonsBy(txdb, by="gene")
saveRDS(genes, file = "exon_by_genes_gencode.v26.GRCh38.rds")


#
#produce counts
#
snps = "datagen/snps"
snpf = list.files(snps, pattern="\\.txt$")
bams = "data/bam"
bamf = list.files(bams, pattern="\\.bam$")

cnts = "data/cnt"
sam = read.table(sprintf("%s/samples.txt", cnts), header=T, as.is=T)

for(i in 1:nrow(sam)){
  sami = grep(sam$ID[i], bamf)
  snpi = grep(sam$ID[i], snpf)
  bam_file = sprintf("%s/%s", bams, bamf[sami]); bam_file
  gene_anno_file = "exon_by_genes_gencode.v26.GRCh38.rds"
  sam_name = sam$ID[i]
  het_snp_file = sprintf("%s/%s", snps, snpf[snpi])
  rcmdi = "Docker/get_TReC_ASReC.R"
  rcmdo = sprintf("Docker/get_TReC_ASReC_%s.Rout", i)
  
  cmd = sprintf("R CMD BATCH '--args %s %s %s %s' %s %s",
  bam_file, gene_anno_file, sam_name, het_snp_file, rcmdi, rcmdo)
  system(cmd)
  message(i)
}

TReC = ASReC_hap1 = ASReC_hap2 = matrix(0, nrow=nrow(gene_info), ncol=nrow(sam))
for(i in 1:nrow(sam)){
  resi = read.table(sprintf("%s.trecase.txt", sam$ID[i]))
  m = match(gnkp, rownames(resi))
  resi = resi[m,]
  write.table(resi, sprintf("%s/%s.trecase.txt", cnts, sam$ID[i]))
  TReC[,i] = resi[,1]
  ASReC_hap1[,i] = resi[,2]
  ASReC_hap2[,i] = resi[,3]
}
write.table(TReC, sprintf("%s/TReC_ex.txt", cnts), row.names=F, col.names=F, quote=F)
write.table(ASReC_hap1, sprintf("%s/ASReC_hap1_ex.txt", cnts), row.names=F, col.names=F, quote=F)
write.table(ASReC_hap2, sprintf("%s/ASReC_hap1_ex.txt", cnts), row.names=F, col.names=F, quote=F)


q("no")


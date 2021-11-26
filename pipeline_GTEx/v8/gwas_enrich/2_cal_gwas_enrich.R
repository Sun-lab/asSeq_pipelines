# just run the simulation example. 

source("gwas_enrich.R")
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)
library(doRNG)

nCore = 11      # number of cores for multi-core computation
registerDoParallel(cores=nCore)

options(mc.cores=nCore)

source_dir = "~/research/GTEx/eQTL_Zhabotynsky_2021/gwas_enrich"
work_dir   = "."

enrich_files = list.files(source_dir)
tissues = gsub(".tsv.gz", "", enrich_files)

gwas = get_GWAS_catalog(work_dir = work_dir)
tab = table(gwas$myPHENO); tab
length(tab)

if(file.exists("gwas_enrich_28_tissues.txt")){
  res_all = fread("gwas_enrich_28_tissues.txt")
}else{
  res_all = NULL
  
  for(i in 1:length(tissues)){
    
    cat("\n------------------------------------------------------\n")
    cat(i, date(), tissues[i])
    cat("\n------------------------------------------------------\n")
    
    dat = fread(file.path(source_dir, enrich_files[i]), data.table=FALSE)
    
    names(dat)[which(names(dat) == "Pos")] = "POS"
    dim(dat)
    dat[1:2,]
    
    dat$Chr = paste0("chr", dat$Chr)
    
    # Run enrichment analyses for all phenotypes
    res0 = run_gwasEnrich_analysis(DATA = dat, work_dir = work_dir,
                                   which_gwas = "gwas_catalog", nBLOCKS = 200, 
                                   verbose = FALSE)
    
    res = foreach(wGWAS = names(tab), .combine = "rbind") %dopar% {
      res1 = run_gwasEnrich_analysis(DATA = dat, work_dir = work_dir,
                                     which_gwas = wGWAS, nBLOCKS = 200,
                                     verbose = FALSE)
      res1
    }
    
    dim(res)
    res[1:2,]
    
    res = rbind(res0, res)
    
    names(res)[2] = "MODEL"
    
    res$GROUP = tissues[i]
    res_all = rbind(res_all, res)
  }
  fwrite(res_all, "gwas_enrich_28_tissues.txt", sep="\t")
}

dim(res_all)
res_all[1:2,]

table(table(res_all$GROUP))
res_all$GROUP = substr(res_all$GROUP, 1, 20)
table(table(res_all$GROUP))
tissues = substr(tissues, 1, 20)
# Plot

res_all$tmp_col = ifelse(res_all$log_enrich_lowJK > 0, 
                         "Significant", "Non-significant")
res_all$wGWAS = as.character(res_all$wGWAS)
lev_wGWAS = sort(unique(res_all$wGWAS)); # lev_wGWAS
lev_wGWAS = c("gwas_catalog",lev_wGWAS[lev_wGWAS != "gwas_catalog"])
lev_wGWAS = rev(lev_wGWAS)
res_all$wGWAS = factor(res_all$wGWAS,levels = lev_wGWAS)

pd = position_dodge(width = 0.75) # control point spread
themes = theme(text = element_text(size = 24),
               panel.background = element_blank(),
               panel.spacing = unit(0.5,"lines"),
               panel.border = element_rect(color = "black",fill = NA,size = 1),
               legend.position = c("none","bottom")[2],
               legend.text = element_text(size = 26))

for(k in 1:7){
  
  res = res_all[res_all$GROUP %in% tissues[(4*(k-1)+1):(4*k)],]
  
  tmp_range = max(abs(res$log_enrich_meanJK)); tmp_range
  tmp_range = 1.2 * c(-1,1) * ifelse(tmp_range >= 1,tmp_range,1); tmp_range
  
  gg = ggplot(data = res,aes(x = wGWAS,y = log_enrich_meanJK,
                             ymin = log_enrich_lowJK,
                             ymax = log_enrich_highJK,
                             group = MODEL)) +
    geom_errorbar(position = pd, size = 1 ,width = 0.1, aes(color = MODEL)) +
    geom_point(position = pd, size = 5,
               aes(shape = tmp_col, color = MODEL, stroke = c(2,1.5)[2])) + 
    geom_point(data = res, position = pd, 
               mapping = aes(x = wGWAS,y = log_enrich, group = MODEL)) +
    facet_grid(~ GROUP) + labs(shape = "Inference", color = "Method") +
    scale_shape_manual(values = c(1,19)) +
    geom_hline(yintercept = 0,linetype = 2) +
    ylab("log(Enrichment)") + xlab("") +
    coord_flip(ylim = tmp_range) + themes + 
    guides(color = guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size = 5)))
  
  fn = file.path(work_dir, sprintf("enrich_trecase_vs_ols_f%d.pdf", k))
  ggsave(filename = fn, plot = gg, width = 15, height = 10, units = "in")
}


gc()

sessionInfo()
q(save="no")

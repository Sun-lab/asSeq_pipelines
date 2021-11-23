# just run the simulation example. 

source("gwas_enrich.R")
library(data.table)
library(ggplot2)

source_dir = "~/research/GTEx/eQTL_Zhabotynsky_2021/gwas_enrich"
work_dir   = "."

enrich_files = list.files(source_dir)
tissues = gsub(".tsv.gz", "", enrich_files)

i = 1
dat = fread(file.path(source_dir, enrich_files[i]), data.table=FALSE)

names(dat)[which(names(dat) == "Pos")] = "POS"
dim(dat)
dat[1:2,]
dat$Chr = paste0("chr", dat$Chr)

# Run enrichment analyses: by phenotype group and across all phenotypes
res = run_gwasEnrich_analysis(DATA = dat, work_dir = work_dir,
                              which_gwas = "gwas_catalog", nBLOCKS = 200, 
                              verbose = FALSE)

gwas = get_GWAS_catalog(work_dir = work_dir)
tab = table(gwas$myPHENO); tab

for(wGWAS in names(tab)){
  cat(sprintf("%s: wGWAS = %s ...\n",date(),wGWAS))
  res1 = run_gwasEnrich_analysis(DATA = dat, work_dir = work_dir,
                                 which_gwas = wGWAS, nBLOCKS = 200,
                                 verbose = FALSE)
  res = rbind(res,res1)
}

dim(res)

fres = rbind(ST_smart_df(res,  MODEL = "CSeQTL"))
#             ST_smart_df(res2, MODEL = "OLS"))

# Plot
pd = position_dodge(width = 0.75) # control point spread
themes = theme(text = element_text(size = 24),
               # axis.text.x = element_text(size = 12),
               # axis.title.y = element_text(size = 20,face = "bold"),
               panel.background = element_blank(),
               panel.spacing = unit(0.5,"lines"),
               panel.border = element_rect(color = "black",fill = NA,size = 1),
               legend.position = c("none","bottom")[2],
               legend.text = element_text(size = 26))


tmp_range = max(abs(fres$log_enrich_meanJK)); tmp_range
tmp_range = 1.2 * c(-1,1) * ifelse(tmp_range >= 1,tmp_range,1); tmp_range
fres$tmp_col = ifelse(fres$log_enrich_lowJK > 0,"Significant","Non-significant")
fres$wGWAS = as.character(fres$wGWAS)
lev_wGWAS = sort(unique(fres$wGWAS)); # lev_wGWAS
lev_wGWAS = c("gwas_catalog",lev_wGWAS[lev_wGWAS != "gwas_catalog"])
lev_wGWAS = rev(lev_wGWAS)
fres$wGWAS = factor(fres$wGWAS,levels = lev_wGWAS)

gg = ggplot(data = fres,aes(x = wGWAS,y = log_enrich_meanJK,
                            ymin = log_enrich_lowJK, 
                            ymax = log_enrich_highJK, group = MODEL)) +
  geom_errorbar(position = pd,size = 1,width = 0.1,aes(color = MODEL)) +
  geom_point(position = pd,size = 5,
             aes(shape = tmp_col,color = MODEL,stroke = c(2,1.5)[2])) + 
  geom_point(data = fres,position = pd,
             mapping = aes(x = wGWAS,y = log_enrich,group = MODEL)) +
  facet_grid(~ GROUP) + labs(shape = "Inference",color = "Cohort") +
  scale_shape_manual(values = c(1,19)) +
  geom_hline(yintercept = 0,linetype = 2) +
  ylab("log(Enrichment)") + xlab("") +
  coord_flip(ylim = tmp_range) +
  themes +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

png_fn = file.path(work_dir,"enrich_test.png")
ggsave(filename = png_fn, plot = gg, width = 15, height = 13, units = "in")
gc()

sessionInfo()
q(save="no")
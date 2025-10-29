setwd('/adpelle1/xqtl-paper-final/main_text/6_AD_xQTL_genes/staging/gene_priorization_table/')

source('gene_prio_utils.R')

# install.packages(c('openxlsx','ggraph'))


out<-'analyses_summary/'

res_adxf<-fread('~/adpelle1/xqtl-resources/data/genes/res_allanalysis_ADloci_overlap_selectedcols.csv.gz')

#all loci (with GWAS pval<1e-5)
unique(res_adxf$gene_name)|>length()#1128 genes
table(unique(res_adxf,by=c('gene_name','Method'))$Method)
unique(res_adxf[Method!='trans_finemapping']$gene_name)|>length()#546 genes

unique(res_adxf$locus_index)|>length()#192 loci 
unique(res_adxf$variant_ID)|>length()#11553 variants
unique(res_adxf[max_variant_inclusion_probability>0.1]$variant_ID)|>length()#697 variants



#genome wide significant only
unique(res_adxf[min_pval<5e-8]$gene_name)|>length()#248 genes genome wide significant
unique(res_adxf[min_pval<5e-8]$locus_index)|>length()#103 loci genome wide significant
unique(res_adxf[min_pval<5e-8]$variant_ID)|>length()#3220 variant genome wide significant
unique(res_adxf[min_pval<5e-8][max_variant_inclusion_probability>0.1]$variant_ID)|>length()#404 variants


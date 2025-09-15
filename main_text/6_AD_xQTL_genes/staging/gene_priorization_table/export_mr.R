setwd('/adpelle1/xqtl-paper/')

source('../../../alexandre-utils/r_utils.R')
source('codes/utilis.R')
source('codes/cb_plot.R')
source('codes/qtl_utils.R',chdir = TRUE)
out<-'analyses_summary/'

res_twas<-fread('analyses_summary/res_AD_XWAS.csv.gz')
table(res_twas$gwas_source)
#for all gene_id
#get file paths
#format: /data/analysis_result/twas/[context].[block].mr_results.tsv.gz
files<-fread('files_twas_folder.txt',col.names = 'file',header = F)
files[,context:=strsplit(file,'\\.')[[1]][1],by='file']
files[,block:=strsplit(file,'\\.')[[1]][2],by='file']
files[,result_type:=strsplit(file,'\\.')[[1]][3],by='file']
fwrite(files,'twas_results_metadata.csv.gz')
files<-fread('twas_results_metadata.csv.gz')

filesf<-files[result_type=='mr_result']
table(filesf$context)
filesfa<-merge(filesf,unique(res_twas[,.(block,gene_ID)]),by = 'block',allow.cartesian = T)
table(filesfa$context)

mr_res<-rbindlist(lapply(unique(filesfa$file),function(f){
  i<-which(unique(filesfa$file)==f)
  if(i%%100==0){
    message(i,'/',length(unique(filesfa$file)))
  }
  genes=filesfa[file==f]$gene_ID
  p<-fp('/data/analysis_result/twas/',f)
  mr<-fread(p)
  if(nrow(mr)>0){
    return(mr[!is.na(meta_pval)][gene_name%in%genes])
    
  }else{
    return(mr)
    
  }
  
  
}),fill=T)

setnames(mr_res,'gene_id','gene_ID')

mr_resf<-mr_res[!is.na(gene_ID)]

#add chr pos of twas
mr_resf<-merge(mr_resf,res_twas[,.(`#chr`,start,end,gene_ID,TADB_start,TADB_end,context)],
      by=c('context','gene_ID'),all.x = T)

#export
fwrite(mr_resf[order(`#chr`,start)],'/adpelle1/export/FunGen_mr.exported.bed.gz')
mr_resf<-fread('/adpelle1/export/FunGen_mr.exported.bed.gz')

table(mr_resf$gwas_source)
#add to mtd
mtd<-fread('all_analysis_summary_tables_metadata.csv')
mtd<-rbind(mtd,data.table(
                          'Data Type'='Gene & GWAS',
                          Cohort='ROSMAP & MSBB & AD',
                          Method='MR',
                          Path='interactive_analysis/adpelle1/export/FunGen_mr.exported.bed.gz'),
           fill=TRUE)
mtd[Method=='MR']
fwrite(mtd,'all_analysis_summary_tables_metadata.csv')


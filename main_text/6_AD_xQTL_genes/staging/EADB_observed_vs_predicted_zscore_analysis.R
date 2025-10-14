setwd('/adpelle1/xqtl-paper-final/')
out<-'main_text/6_AD_xQTL_genes/staging/APOElocus/'
dir.create(out)
source('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/gene_prio_utils.R')

library(pecotmr)

#get the xQTL in APOE region
resapoe<-fread(fp(out,'xqtl_in_APOE50kbregion.csv.gz')) #generated in merge_APOE_cs.R

#remove all
MicQTL_unioncs_variants<-resapoe[uni.locus_id%in%c('loc_1608','loc_1610')]$variant_ID|>unique()|>str_remove('chr')
saveRDS(MicQTL_unioncs_variants,fp(out,'MicQTL_cs_union_3_and_4_variants.rds'))
# this have been generated taking union of all CS from finemapping. i.e. if 2 CS of different dataset overlap, we take the union of them as a single CS


#get slalom QC imputed gwas z scores####
#/!\ here I use slalom instead of rss_qc, and made some specific adaptation to account for the huge zscores for this region leading to instability
##if zscore =Inf, I transformed as the max value
#- in slalom, I used r2_threshold =0.2 instead of default 0.6 for outlier detection because otherwise APOE4 and others important variant, the observed zscore are removed and then imputed
#rationale being it is important to conserve the observed zscore for those variants if we want to be more accurate in imputation
#we use as gwas summ stat EADB only, because coming from a single european cohort (no metaanalysis) of AD and Control case, ensuring more homogenity and alignment with our LD panel

impute_opts=list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)
qc_method='slalom'
gwas='AD_Bellenguez_EADB_2022'

gwasmtd<-fread('/data/GWAS/GWAS_sumstat_meta_Apr_9.tsv')


ld_meta<-fread('/data/resource/ADSP_R4_EUR_LD/ld_meta_file_apoe.tsv')
ld_meta#should be updated because path not correct (need to rm 'chr19/')
ldblock<-load_LD_matrix('../xqtl-paper/resources/ld_meta_file_apoe.tsv',
                        region = data.frame(chr='chr19',
                                            start=42346101,
                                            end=46842901))




sumstat_path='/data/GWAS/AD_GWAS/EADB_core.tsv.gz'
column_file_path='/data/GWAS/column_mapping_file/EADB.yml'

#get summstats
rss_input <- load_rss_data(
  sumstat_path = sumstat_path, column_file_path = column_file_path
)
rss_input$n

# Preprocess the input data
preprocess_results <- rss_basic_qc(rss_input$sumstats, ldblock)

#if zscore =Inf, transform as the max value
max(preprocess_results$sumstats$z[!is.infinite(preprocess_results$sumstats$z)])
if(any(is.infinite(preprocess_results$sumstats$z))){
  message(sum(is.infinite(preprocess_results$sumstats$z)),' variants have infinite zscore, exclude them')
  preprocess_results$sumstats<-preprocess_results$sumstats[!is.infinite(preprocess_results$sumstats$z),]
  preprocess_results$LD_mat=preprocess_results$LD_mat[preprocess_results$sumstats$variant_id,preprocess_results$sumstats$variant_id]
}

# Perform Slalom quality control

# qc_results <- summary_stats_qc(preprocess_results$sumstats, list(combined_LD_matrix=preprocess_results$LD_mat), n = rss_input$n,
#                                var_y = rss_input$var_y,
#                                method = qc_method) #cannot change the r2_threshold of slalom here so used directly the slalom funciton

res_slalom <- slalom(zScore = preprocess_results$sumstats$z,
                     LD_mat = preprocess_results$LD_mat,r2_threshold=0.2)
qc_results<-list()
qc_results$sumstats<-merge(data.table(preprocess_results$sumstats),
                           data.table(res_slalom$data,keep.rownames = 'variant_id')[,.(variant_id,nlog10p_dentist_s)],
                           by='variant_id')
qc_results$LD_mat<-preprocess_results$LD_mat[qc_results$sumstats$variant_id,qc_results$sumstats$variant_id]
qc_results$outlier_number<-nrow(data.table(res_slalom$data,keep.rownames = 'variant_id')[(outliers)])


# Perform imputation
impute_results <- raiss(ldblock$ref_panel, qc_results$sumstats, ldblock$combined_LD_matrix,
                        rcond = impute_opts$rcond,
                        R2_threshold = impute_opts$R2_threshold,
                        minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb)

summstats<-data.table(impute_results$result_filter)

#if imputed zscore!=orignial and original zscore>15, keep original
commonzs<-merge(
  data.table(preprocess_results$sumstats)[,.(variant_id,z)],
  summstats[,.(variant_id,z)], 
  by='variant_id')

var_toswitch<-commonzs[z.x>z.y&z.x>15]$variant_id
message(length(var_toswitch),' variants to recover observed zscore')
if(length(var_toswitch)>0){
  summstats[var_toswitch,z:=data.table(preprocess_results$sumstats)[var_toswitch,on='variant_id']$z,on='variant_id']
  
}


#add variants annots
#add if from cs3, cs4 and leading variants based on pip
resapoe[,top_snp:=PIP==max(PIP),by='uni.locus_id']
resapoe[uni.locus_id=='loc_1608'&top_snp]
resapoe[,variant_id2:=str_remove(variant_ID,'chr')]
summstats[,variant_anno:=ifelse(variant_id=='19:44908684:T:C','APOE4',
                               ifelse(variant_id=='19:44908822:C:T','APOE2',
                                      ifelse(variant_id==resapoe[uni.locus_id=='loc_1608'&top_snp]$variant_id2,'cs3 lead',
                                             ifelse(variant_id%in%resapoe[uni.locus_id=='loc_1608']$variant_id2,'cs3',
                                                    ifelse(variant_id==resapoe[uni.locus_id=='loc_1610'&top_snp]$variant_id2,'cs4 lead',
                                                    ifelse(variant_id%in%resapoe[uni.locus_id=='loc_1610']$variant_id2,'cs4',''))))))]

#annot with correl with cs3 lead and apoe4
summstats[,apoe4_r:=ldblock$combined_LD_matrix[variant_id,
                                  '19:44908684:T:C']]
summstats[,apoe2_r:=ldblock$combined_LD_matrix[variant_id,
                                  '19:44908822:C:T']]

summstats[,cs3_lead_r:=ldblock$combined_LD_matrix[variant_id,
                                                  resapoe[uni.locus_id%in%c('loc_1608')&top_snp]$variant_id2[1]]]

fwrite(summstats,file.path(out,'eadb_zscore_slalom_qc_APOEblock.csv.gz'))
summstats<-fread(file.path(out,'eadb_zscore_slalom_qc_APOEblock.csv.gz'))

#predicted zscore vs observed for the variants of interest####
#we try both by imputing based on variants at +/-50kb of APOE, and then the full block
# 1) impute only based on the +/-50kb APOE region ####

apoereg<-get_gene_info(gene_name = 'APOE')
start.reg=apoereg$gene_info$gene_start-50000
end.reg=apoereg$gene_info$gene_end+50000
reg=paste('chr19',start.reg,end.reg,sep='-')


summstatsf<-summstats[pos>start.reg&pos<end.reg]

variants_to_include<-unique(summstatsf[order(pos)]$variant_id)


#imputation of variants of interest,=> excluded them from the observed zscore
summstatsff<-summstatsf[!variant_id%in%MicQTL_unioncs_variants]



  zscores=summstatsff[,.(chrom=seqid(variant_id),
                         pos=pos(variant_id),
                         variant_id,
                         z=z)]|>as.data.frame()
  
  zscores$A2<-data.table(ldblock$ref_panel)[zscores$variant_id,on='variant_id']$A2
  zscores$A1<-data.table(ldblock$ref_panel)[zscores$variant_id,on='variant_id']$A1
  
  imputed_zscores<-raiss(ref_panel = data.table(ldblock$ref_panel)[variants_to_include,on='variant_id']|>as.data.frame(),
                         known_zscores =zscores,
                         LD_matrix=ldblock$combined_LD_matrix[variants_to_include,variants_to_include],
                         R2_threshold=0.4,minimum_ld = 0.5)
  
  res_imput<-data.table(imputed_zscores$result_nofilter)[,pass_qc:=variant_id%in%imputed_zscores$result_filter$variant_id]
 #merge imput and predicted
   setnames(res_imput,'z','predicted.zscore')
   setnames(summstatsf,'z','observed.zscore')
   
  
  imputobs<-merge(res_imput,summstatsf,by=c('variant_id'))
  imputobs[predicted.zscore!=observed.zscore]|>nrow()#61
  imputobs[variant_id%in%MicQTL_unioncs_variants]
#Zoom on +/-3  zscore of the cs3 leading variant
imputobs[,zoom:=predicted.zscore<max(predicted.zscore[variant_anno=='cs3 lead'])+3&predicted.zscore>min(predicted.zscore[variant_anno=='cs3 lead'])-3]

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=variant_anno))+theme_bw()

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=abs(cs3_lead_r)))+theme_bw()

ggplot(imputobs[(zoom)][])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=abs(apoe4_r)))+theme_bw()+
  scale_color_gradient2(midpoint = 0.1,low = 'blue',high = 'red',mid = 'grey')

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=apoe2_r))+theme_bw()
  
#Zoom  on +/- 3 zscore of all variants of interest
imputobs[,zoom:=predicted.zscore<max(predicted.zscore[variant_id%in%MicQTL_unioncs_variants])+3&predicted.zscore>min(predicted.zscore[variant_id%in%MicQTL_unioncs_variants])-3]

ggplot(imputobs[(zoom)])+
geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=variant_anno))+theme_bw()

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=cs3_lead_r))+theme_bw()

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=apoe4_r))+theme_bw()

ggplot(imputobs[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=apoe2_r))+theme_bw()


#2) imputing with all block####

variants_to_include<-unique(summstats[order(pos)]$variant_id)


#imputation of variants of interest,=> excluded them from the observed zscore
summstatsff<-summstats[!variant_id%in%MicQTL_unioncs_variants]



zscores=summstatsff[,.(chrom=seqid(variant_id),
                       pos=pos(variant_id),
                       variant_id,
                       z=z)]|>as.data.frame()

zscores$A2<-data.table(ldblock$ref_panel)[zscores$variant_id,on='variant_id']$A2
zscores$A1<-data.table(ldblock$ref_panel)[zscores$variant_id,on='variant_id']$A1

imputed_zscores<-raiss(ref_panel = data.table(ldblock$ref_panel)[variants_to_include,on='variant_id']|>as.data.frame(),
                       known_zscores =zscores,
                       LD_matrix=ldblock$combined_LD_matrix[variants_to_include,variants_to_include],
                       R2_threshold=0.4,minimum_ld = 0.5)

res_imput<-data.table(imputed_zscores$result_nofilter)[,pass_qc:=variant_id%in%imputed_zscores$result_filter$variant_id]
#merge imput and predicted
setnames(res_imput,'z','predicted.zscore')
setnames(summstats,'z','observed.zscore')


imputobsfull<-merge(res_imput,summstats,by=c('variant_id'))
imputobsfull[predicted.zscore!=observed.zscore]|>nrow()#61
imputobsfull[variant_id%in%MicQTL_unioncs_variants]
#Zoom on +/-3  zscore of the cs3 leading variant
imputobsfull[,zoom:=predicted.zscore<max(predicted.zscore[variant_id%in%MicQTL_unioncs_variants])+3&predicted.zscore>min(predicted.zscore[variant_id%in%MicQTL_unioncs_variants])-3]

ggplot(imputobsfull[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=variant_anno))+theme_bw()

ggplot(imputobsfull[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=cs3_lead_r))+theme_bw()

ggplot(imputobsfull[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=apoe4_r))+theme_bw()

ggplot(imputobsfull[(zoom)])+
  geom_point(aes(x=abs(predicted.zscore),y=abs(observed.zscore),col=apoe2_r))+theme_bw()


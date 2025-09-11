
#Gene locus summary
setwd('/adpelle1/xqtl-paper/')
install.packages(c('openxlsx','ggraph'))

source('../../../alexandre-utils/r_utils.R')
source('codes/qtl_utils.R')

out<-'analyses_summary/'
dir.create(out)

#0)prep integration####
#harmonization of contexts names
# contexts<-fread('contexts_metadata.csv')
#merge with coloc context
# coloccont<-fread('colocname_to_contexts.tsv')
# coloccont[context=='DLPFC_Klein_gpQTL',context:='DLPFC_Klein_gpQTL_unadjusted']
# fwrite(coloccont,'colocname_to_contexts.tsv')
# contexts<-merge(contexts,coloccont,all.x = T,by='context')
# setdiff(coloccont$context,contexts$context)
# fwrite(contexts,'contexts_metadata.csv')
# #merge with Qr context name
# #on excel

#I) Prep/harmionized the exported tables for each methods/resource####
#for each, long the table at variant context, gene level and create an unique locus-event identifier

mtd<-fread('all_analysis_summary_tables_metadata.csv',header = T)
contexts<-fread('contexts_metadata.csv')
# setnames(contexts,'context_hi','context_broad')
# mtd<-merge(mtd,unique(contexts),all.x = T,by='context')
# fwrite(mtd[,.(`Data Type`,Cohort,Modality,context,context_broad,Method,Path)],'all_analysis_summary_tables_metadata.csv')
table(mtd$Method)
# AD_GWAS_finemapping           APOE interaction                      Coloc                 ColocBoost                      ctwas 
# 8                         23                          2                          3                          1 
# fSuSiE_finemapping                         LR           msex interaction  multi_context_finemapping     multi_gene_finemapping 
# 10                         35                         35                          2                         21 
# QR single_context_finemapping          trans_finemapping                       twas 
# 22                         48                         36                          1 

# 1) the finemapping tables####
#  Finemapping: singlecontext, multicontexts, fsusie, multi gene, transqtl
# - singlecontexts####
res_sc<-rbindlist(lapply(file.path('/data',mtd[Method=='single_context_finemapping']$Path),function(f)fread(f)),fill = T)
res_sc[!is.na(cat)]
res_sc[,context:=str_extract(event_ID[1], '^.+?(?:_chr|_ENSG|_gp_)')%>%
         gsub('(_chr|_ENSG|_gp_)', '', .),by='event_ID']
table(res_sc$context)
#create locus_context_id for cs95>cs70>cs50
res_sc[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_scf<-res_sc[(hit)]
res_scf[,
     credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                         ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                paste0('cs50_',cs_coverage_0.5)))]
res_scf[,locuscontext_id:=paste(event_ID,credibleset,sep='_')]
res_scf[,region:=gene_ID]
fwrite(res_scf,fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz'))
res_scf<-fread(fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz'))
unique(res_scf,by='variant_ID')[alt(variant_ID)==a2]
table(res_scf$context)

mtd[Method=='single_context_finemapping',summary_file:=fp(out,'res_all_single_context_finemapping_cs50orgreater.csv.gz')]

# - multicontexts####
res_mc<-rbindlist(lapply(file.path('/data',mtd[Method=='multi_context_finemapping']$Path),function(f)fread(f)),fill = T)
res_mc[lfsr=='']
#first rm thos without lfsr
res_mcf<-res_mc[!lfsr=='']

#one row per context
multi<-res_mcf[,.(context=strsplit(event_ID,';')[[1]],
                                         conditional_effect=strsplit(conditional_effect,';')[[1]]|>as.numeric(),
                                         lfsr=strsplit(lfsr,';')[[1]]|>as.numeric()),
               by=c('event_ID','gene_ID')]|>unique()
res_mcf<-merge(res_mcf[,-c('conditional_effect','lfsr')],multi,by=c('event_ID','gene_ID'),allow.cartesian = T)

#create locus_context_id for cs95>cs70>cs50
res_mcf[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_mcf<-res_mcf[(hit)]
res_mcf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_mcf[,locuscontext_id:=paste(event_ID,gene_ID,credibleset,sep='_')]
res_mcf[,region:=gene_ID]
max(res_mcf$lfsr)
fwrite(res_mcf,fp(out,'res_all_multi_context_finemapping_cs50orgreater_lfsr0.01.csv.gz'))

mtd[Method=='multi_context_finemapping',summary_file:=fp(out,'res_all_multi_context_finemapping_cs50orgreater_lfsr0.01.csv.gz')]

# - fsusie##### 
file.exists(file.path('/data',mtd[Method=='fSuSiE_finemapping']$Path))
res_fs<-rbindlist(lapply(file.path('/data',mtd[Method=='fSuSiE_finemapping']$Path),function(f)fread(f)),fill = T)
res_fs[,1:17]
colnames(res_fs)
res_fs[variant_ID=='chr16:30067171:G:C']$cs_id
#create locus_context_id for cs95>cs70>cs50
res_fs[,hit:=cs_coverage_0.95!=0]
res_fs[(hit)]#everything is cs95
res_fsf<-res_fs
# res_fsf<-res_fs[(hit)]
# res_fsf[,
#         credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
#                             ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
#                                    paste0('cs50_',cs_coverage_0.5)))]
res_fsf[,locuscontext_id:=cs_id]
res_fsf[,region:=region_ID]
res_fsf[,context:=event_ID]
res_fsf[,chr:=seqid(variant_ID)]
res_fsf[,pos:=pos(variant_ID)]

setdiff(res_fsf$context,mtd$context)
res_fsf<-res_fsf[,-'grid_position']
colnames(res_fsf)
res_fsf[,top_effect_coord:=strsplit(epi_mark_effects[1],split = ';')[[1]]|>as.numeric()|>abs()|>which.max(),by='locuscontext_id']
res_fsf[,top_effect:=strsplit(epi_mark_effects[1],split = ';')[[1]][top_effect_coord]|>as.numeric(),by='locuscontext_id']
res_fsf[,n_epi_marks:=strsplit(epi_mark_effects[1],split = ';')[[1]]|>as.numeric()|>length(),by='locuscontext_id']
res_fsf[,-'grid_positions']
summary(res_fsf$n_epi_marks)
summary(res_fsf[str_detect(context,'ATAC')]$n_epi_marks)
table(unique(res_fsf,by='locuscontext_id')$context)
res_fsf[str_detect(context,'ATAC')][,.(epi_mark_names)]
fwrite(res_fsf,fp(out,'res_all_fsusie_finemapping_cs95.csv.gz'))
res_fsf<-fread(fp(out,'res_all_fsusie_finemapping_cs95.csv.gz'))
mtd[Method=='fSuSiE_finemapping',summary_file:=fp(out,'res_all_fsusie_finemapping_cs95.csv.gz')]

# res_fsf<-UnifyLoci(res_fsf,variant_col = 'variant_ID',locus_col = 'locuscontext_id',group.by = 'region_ID')

# - multi gene####
res_mg<-rbindlist(lapply(file.path('/data',mtd[Method=='multi_gene_finemapping']$Path),function(f)fread(f)[,Path:=f]),fill = T)
res_mg[lfsr=='']
#first rm thos without lfsr
res_mgf<-res_mg[!lfsr=='']
res_mgf[,region:=gene_ID]
#one row per gene
multi<-res_mgf[,.(gene_ID=strsplit(event_ID,';')[[1]]|>str_extract('ENSG[0-9]+'),
                  conditional_effect=strsplit(conditional_effect,';')[[1]]|>as.numeric(),
                  lfsr=strsplit(lfsr,';')[[1]]|>as.numeric()),
               by=c('event_ID','region')]|>unique()
res_mgf<-merge(res_mgf[,-c('conditional_effect','lfsr','gene_ID')],multi,by=c('event_ID','region'),allow.cartesian = T)

#create locus_context_id for cs95>cs70>cs50
res_mgf[,hit:=cs_coverage_0.95!=0]

res_mgf<-res_mgf[(hit)]
res_mgf[,
        credibleset:=paste0('cs95_',cs_coverage_0.95)]
#need add context
res_mgf[,context:=strsplit(basename(Path[1]),'\\.')[[1]][1],by='Path']
setdiff(res_mgf$context,contexts$context) 

res_mgf[str_detect(context,'Oli|Exc|Inh|Ast|Mic|OPC|monocyte|DLPFC|AC|PCC'),context:=str_remove(context,'^ROSMAP_')]
res_mgf[str_detect(context,'MSBB'),context:=paste0('BM_',str_extract(context,'[0-9]+'),'_MSBB_eQTL')]
setdiff(res_mgf$context,contexts$context) 

res_mgf[,locuscontext_id:=paste(context,'multigene',event_ID,credibleset,sep='_')]

fwrite(res_mgf,fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz'))
res_mgf<-fread(fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz'))

mtd[Method=='multi_gene_finemapping',summary_file:=fp(out,'res_all_multi_gene_finemapping_cs50orgreater.csv.gz')]


# - the transQTL####
res_ts<-rbindlist(lapply(file.path('/data',mtd[Method=='trans_finemapping']$Path),function(f)fread(f)),fill = T)
res_ts[!is.na(cat)]
res_ts[,context:=str_extract(event_ID[1], '^.+?(?:_chr|_ENSG|_gp_)')%>%
         gsub('(_chr|_ENSG|_gp_)', '', .),by='event_ID']
table(res_ts$context)
#create locus_context_id for cs95>cs70>cs50
res_ts[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_tsf<-res_ts[(hit)]
res_tsf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_tsf[,locuscontext_id:=paste(event_ID,gene_ID,credibleset,sep='_')]
res_tsf[,region:=gene_ID]
res_tsf[,gene_ID:=str_extract(event_ID,'ENSG[0-9]+')]

fwrite(res_tsf,fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz'))
res_tsf<-fread(fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz'))

mtd[Method=='trans_finemapping',summary_file:=fp(out,'res_all_transgene_single_context_finemapping_cs50orgreater.csv.gz')]

# the gwas results: consolidate gwas locus with min corr 0.5 ####
res_gw<-rbindlist(lapply(file.path('/data',mtd[Method=='AD_GWAS_finemapping']$Path),function(f)fread(f)),fill = T)
res_gw[!is.na(cat)]
res_gw[,context:=event_ID]
res_gw[,region:=gene_ID]
res_gw[,study:=event_ID]
#
#consolidate gwas locus with min corr 0.5 res_gwf<-res_gw[(hit)]
res_gw[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]
res_gwf<-res_gw[(hit)]
fwrite(res_gwf,fp(out,'res_all_single_gwas_finemapping_cs50orgreater.csv.gz'))

# 2) The non linear tables####
# the QR table####
#two pvalue: bonferonni wich is adjusted for each genes, while qvalue is adjusted by number of variants tested within a gene
#for 1
res<-fread("/data/analysis_result/marginal_significant_qtl/cis_association/KNIGHT/eQTL/Brain/QR/Knight_quantile_eQTL.cis_pairs.significant_qtl.filtered_bonferroni_BH_adjusted.tsv.gz")
res[,.(molecular_trait_object_id,variant_id,pvalue,qvalue,p_bonferroni_adj,coef_heter)][coef_heter>-0.6931472]$molecular_trait_object_id|>unique()|>length()
res[,.(molecular_trait_object_id,variant_id,pvalue,qvalue,p_bonferroni_adj,coef_heter)]$molecular_trait_object_id|>unique()|>length()


#need add LR/QR correlation test to separate real non linear effect from linear
#Pending Anjing update
#update in analysis_result/quantile_qtl/export/pure_qr.qtl.bonf_fdr_add_xi.classification.tsv.gz
#edit mtd
mtd[Method=='QR',Path:='analysis_result/quantile_qtl/export/pure_qr.qtl.bonf_fdr_add_xi.classification.tsv.gz']
fwrite(mtd,'metadata_analysis.csv')

res_qr<-rbindlist(lapply(file.path('/data',unique(mtd[Method=='QR']$Path)),function(f){
  
 fread(f)

  
}))
res_qr<-res_qr[p_bonferroni_adj<0.05][,.(chrom,pos,ref,alt,molecular_trait_object_id,variant_id,maf,cis_start,cis_end,feature_tss,feature_tes,pvalue ,qvalue,coef_heter,context,n_variants,p_bonferroni_adj,xi,xi_pval,classification)]

res_qr[,gene_ID:=str_extract(molecular_trait_object_id,'ENSG[0-9]+')]
res_qr[,variant_ID:=str_replace_all(variant_id,'_',':')]
setnames(res_qr,'context','context_qr')
res_qr<-merge(res_qr,contexts[,.(context,context_qr,context_hi)],by = 'context_qr',all.x = T)
res_qr[is.na(context)]
#n genes
unique(res_qr$gene_ID)

fwrite(res_qr,fp(out,'res_qr_summ.csv.gz'))
res_qr<-fread(fp(out,'res_qr_summ.csv.gz'))

mtd[Method=='QR',summary_file:=fp(out,'res_qr_summ.csv.gz')]
mtd<-unique(mtd,by=c('Path'))
mtd[Method=='QR',Cohort:='ROSMAP, MSBB, KNIGHT']
mtd[Method=='QR',`Data Type`:='eQTL & pQTL']
mtd[Method=='QR',`Modality`:='-']

fwrite(mtd,'metadata_analysis.csv')

# the interactions tables####
#have been corrected 
#for 1
interactions<-fread('/data/analysis_result/marginal_significant_qtl/cis_association/KNIGHT/eQTL/Brain/interaction/msex/xqtl_protocol_data.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.bed.bed.cis_pairs.significant_qtl.q_bonferroni_min_adjusted_events_qvalue.tsv.gz')
max(interactions$qvalue_interaction)
#for all
#msex
msex_summ<-rbindlist(lapply(file.path('/data',mtd[Method=='msex interaction']$Path),function(f){
  file=list.files(f,pattern='q_bonferroni_min_adjusted_events_qvalue.tsv.gz',full.names = T)
  if(length(file)==1){
    return(fread(file)[,Path:=str_remove(f,'/data/')])
    
  }else{
    return(data.table(Path=str_remove(f,'/data/')))
  }

  
}),fill = T)
msex_summ[,variant_id:=str_replace_all(variant_id,'_',':')]
msex_summ[,gene_ID:=str_extract(molecular_trait_id,'ENSG[0-9]+')]
msex_summ[,variant_ID:=str_replace_all(variant_id,'_',':')]
msex_summ<-merge(msex_summ,mtd[,.(Path,context,context_hi)],by = 'Path')
msex_summ<-msex_summ[!is.na(variant_id)]
fwrite(msex_summ,fp(out,'res_msex_interaction_summ.csv.gz'))

mtd[Method=='msex interaction',summary_file:=fp(out,'res_msex_interaction_summ.csv.gz')]


#define locus: lead SNP and variants signif in its block with r>0.1
#NOT RUN for now, if needed later
#for one block
# block<-'chr17_42087601_45383525'
# ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
#                         region = data.frame(chr=seqid(block),
#                                             start=start(block),
#                                             end=end(block)))$combined_LD_matrix
# msex_summf<-msex_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
# 
# table(msex_summf$context)
# loci<-getLocusID(msex_summf[context=='DLPFC_Bennett_pQTL'],ldblock,order.by = 'pvalue_msex_interaction',group.by = 'gene_ID')
# table(loci$locus_id)
# 
# 
# #for all
# library(parallel)
# msex_summ<-msex_summ[!is.na(variant_id)]
# ldmeta<-fread('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv')
# blocks<-bed_inter(unique(msex_summ[,.(seqid(variant_id),pos(variant_id)-1,pos(variant_id),variant_id)]),ldmeta,select = 8)[[1]]|>unique()
# ldmetaf<-ldmeta[path%in%blocks]
# ldmetaf[,block:=paste(`#chrom`,start,end,sep='_')]
# msex_summ<-rbindlist(mclapply(unique(ldmetaf$block),function(block){
#   message(block)
#   ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
#                           region = data.frame(chr=seqid(block),
#                                               start=start(block),
#                                               end=end(block)))$combined_LD_matrix
#   summf<-msex_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
#   
#   if(nrow(summf)>0){
#     summf<-getLocusID(summf,ldblock,order.by = 'pvalue_msex_interaction',group.by = c('gene_ID','context'))
#     
#   }else{
#     summf<-data.table()
#   }
# 
#   return(summf)
# },mc.cores = 1),fill = T)

#	APOE interaction
apoe_summ<-rbindlist(lapply(file.path('/data',mtd[Method=='APOE interaction']$Path),function(f){
  file=list.files(f,pattern='q_bonferroni_min_adjusted_events_qvalue.tsv.gz',full.names = T)
  if(length(file)==1){
    return(fread(file)[,Path:=str_remove(f,'/data/')])
    
  }else{
    return(data.table(Path=str_remove(f,'/data/')))
  }
  
  
}),fill = T)
apoe_summ
apoe_summ<-merge(apoe_summ,mtd[,.(Path,context,context_hi)],by = 'Path')
apoe_summ[,gene_ID:=str_extract(molecular_trait_id,'ENSG[0-9]+')]
apoe_summ[,variant_ID:=str_replace_all(variant_id,'_',':')]

apoe_summ<-apoe_summ[!is.na(variant_id)]
fwrite(apoe_summ,fp(out,'res_APOE_interaction_summ.csv.gz'))
apoe_summ<-fread(fp(out,'res_APOE_interaction_summ.csv.gz'))
apoe_summ[pvalue_APOE_interaction]
mtd[Method=='APOE interaction',summary_file:=fp(out,'res_APOE_interaction_summ.csv.gz')]

#define locus: lead SNP and variants signif in its block with r>0.1
#NOT RUN for now, if needed later
# ldmeta<-fread('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv')
# blocks<-bed_inter(unique(apoe_summ[,.(seqid(variant_id),pos(variant_id)-1,pos(variant_id),variant_id)]),ldmeta,select = 8)[[1]]|>unique()
# ldmetaf<-ldmeta[path%in%blocks]
# ldmetaf[,block:=paste(`#chrom`,start,end,sep='_')]
# apoe_summ<-rbindlist(mclapply(unique(ldmetaf$block),function(block){
#   message(block)
#   ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
#                           region = data.frame(chr=seqid(block),
#                                               start=start(block),
#                                               end=end(block)))$combined_LD_matrix
#   summf<-apoe_summ[str_remove(variant_id,'chr')%in%rownames(ldblock)]
#   if(nrow(summf)>0){
#     summf<-getLocusID(summf,ldblock,order.by = 'pvalue_msex_interaction',group.by = c('gene_ID','context'))
#     
#   }else{
#     summf<-data.table()
#   }
#   
#   return(summf)
# },mc.cores = 2))


# 3)the colocs####
table(mtd$Method)
# - xQTLxADs####
res_c<-rbindlist(lapply(file.path('/data',mtd[Method=='ColocBoost'][Modality=='AD_xQTL_colocalization']$Path),
                        function(f)fread(f)),fill=TRUE)
res_c[is.na(gene_ID),gene_ID:=region_ID]
res_c[,gwas_source:=str_extract(event_ID,'AD_[A-Za-z0-9_]+$')]
res_c$gwas_source|>table()
res_c[gwas_source=='AD_Bellenguez',gwas_source:='AD_Bellenguez_2022']
#transform in long format with column gene_id, study, variant_id, vcp,
res_c2<-merge(res_c[,.(
                        variant_ID,
                        vcp,
                        cos_npc,
                        min_npc_outcome),
                      by=.(gene_ID,cos_ID,gwas_source)],
                    res_c[,.(
                      event_ID=strsplit(event_ID,'\\; ')[[1]],
                      coef=as.numeric(strsplit(coef,'\\;')[[1]])),
                       by=.(gene_ID,cos_ID,gwas_source)],
              allow.cartesian=TRUE)
res_c2[event_ID=='AD_Bellenguez',event_ID:='AD_Bellenguez_2022']

res_c2[,top_variant:=vcp==max(vcp),by=.(gene_ID,cos_ID,gwas_source)]
unique(res_c2[(top_variant)],by=c('cos_ID','gene_ID'))
#res_c<-merge(res_c,res_cRu[,.(cos_id,gene_id=gene,top_variable,purity)],by=c('gene_id','cos_id'))
res_c2[,locuscontext_id:=paste(cos_ID,gwas_source,gene_ID,sep='_')]

res_c2[,n.variant:=length(unique(variant_ID)),by=.(locuscontext_id)]
#add the matching context name with fp
res_c2[,context_coloc:=str_remove(event_ID,gene_ID[1])|>str_remove('(_|:)$')|>str_remove('adjusted_gp_[0-9]+|P[0-9]+')|>str_remove('_\\|[A-Z0-9]+$')|>str_remove('chr[0-9]+__[A-Z0-9]+')|>str_remove('_chr[0-9]+:[0-9]+:[0-9]+:clu_[0-9]+_[+-]:[A-Z]+')|>str_remove('chr[0-9]+__')|>str_remove('_\\|')|>str_remove('_$'),by='gene_ID']
unique(res_c2$context_coloc)|>sort()|>cat(sep = '\n')
setdiff(res_c2$context_coloc,contexts$context_coloc)#OK

res_c<-merge(res_c2[,-c('context')],unique(contexts[context_coloc!=''&!str_detect(context,'_u_|_p_')][,.(context_coloc,context)]),by='context_coloc',all.x = T)
res_c[is.na(context)]#ok
table(res_c[is.na(context)]$context_coloc)


fwrite(res_c,fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))


# - ADs #####
res_cad<-fread(file.path('/data',mtd[Method=='ColocBoost'][Modality=='AD_meta_colocalization']$Path))
 #transform in long format with column gene_id, study, variant_id, vcp,
res_cad2<-merge(res_cad[,.(
                              variant_ID,
                              vcp,
                              cos_npc,
                              min_npc_outcome
                            ),
                            by=.(cos_ID)],
                res_cad[,.(
                  event_ID=strsplit(event_ID,'\\; ')[[1]],
                          coef=as.numeric(strsplit(coef,'\\;')[[1]])),
                          
                          by=.(cos_ID)],
                allow.cartesian=TRUE)


res_cad2[event_ID=='AD_Bellenguez',event_ID:='AD_Bellenguez_2022']

res_cad2[event_ID=='AD_Bellenguez_EADB',event_ID:='AD_Bellenguez_EADB_2022']
res_cad2[,context:=event_ID]
res_cad2[,gwas_source:=event_ID]

res_cad2[,top_variant:=vcp==max(vcp),by=.(cos_ID,gwas_source)]
res_cad2[,locuscontext_id:=cos_ID]


fwrite(res_cad2,fp(out,'res_coloc_meta_AD.csv.gz'))

# - fsusie Coloc####
res_cfs<-rbindlist(lapply(file.path('/data',mtd[Method=='Coloc']$Path),function(f)fread(f)),fill=TRUE,use.names = T)

res_cfs[event_ID=='ROSMAP_mQTL',event_ID:='ROSMAP_DLPFC_mQTL']
res_cfs[event_ID=='ROSMAP_haQTL',event_ID:='ROSMAP_DLPFC_haQTL']
res_cfs[,context:=event_ID]
unique(res_cfs$context)

unique(res_cfs$CoS)|>length()
res_cfs[is.na(CoS)]
setdiff(res_cfs$context,contexts$context)

res_cfs[,gwas_source:=AD]
res_cfs$gwas_source|>table()
#keep only those in CoS
res_cfsf<-res_cfs[!is.na(CoS)]
setnames(res_cfsf,'CoS','cos_ID')
#setnames(res_cfsf,'variant_id','variant_ID')

res_cfsf[,top_variant:=SNP_PPH4==max(SNP_PPH4),by=.(cos_ID,gwas_source)]
unique(res_cfsf[(top_variant)],by=c('cos_ID','gwas_source'))
#res_c<-merge(res_c,res_cRu[,.(cos_id,gene_id=gene,top_variable,purity)],by=c('gene_id','cos_id'))
res_cfsf[,locuscontext_id:=paste(TADB_region,gwas_source,context,coloc_index,sep='_')]

res_cfsf[,n.variant:=length(unique(variant_ID)),by=.(locuscontext_id)]
res_cfsf[!str_detect(variant_ID,'^chr'),variant_ID:=paste0('chr',variant_ID)]

res_cfsf[,chr:=seqid(variant_ID)]

fwrite(res_cfsf,fp(out,'res_coloc_AD_epiQTL.csv.gz'))

#  1-3) AD locis Extension/harmonization with FP and ADmeta coloc ####
#run gwas_cs_extension.R

#Finemapping
res_gwf<-fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#create locus_context_id for cs95>cs70>cs50

res_gwf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_gwf[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]

#unified locus
res_gwf[,chr:=seqid(variant_ID)]
res_gwf<-UnifyLoci(res_gwf,variant_col = 'variant_ID',
                   locus_col = 'locuscontext_id',group.by = 'chr')
unique(res_gwf$uni.locus_id)#186
unique(res_gwf[str_detect(credibleset,'cs95')]$uni.locus_id)|>length()#92

res_gwf[,context:=event_ID]
res_gwf[,gwas_source:=event_ID]

fwrite(res_gwf,fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz'))
res_gwf<-fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz'))

mtd[Method=='AD_GWAS_finemapping',summary_file:=fp(out,'res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz')]

# #QC: how many loci without extension?
# res_gwfo<-res_gwf[!is.na(cs_coverage_0.95_original)][cs_coverage_0.5_original!=0|cs_coverage_0.7_original!=0|cs_coverage_0.95_original!=0][,-c('uni.locus_id')]
# res_gwfo[,
#         credibleset:=ifelse(cs_coverage_0.95_original>0,paste0('cs95_',cs_coverage_0.95_original),
#                             ifelse(cs_coverage_0.7_original>0,paste0('cs70_',cs_coverage_0.7_original),
#                                    paste0('cs50_',cs_coverage_0.5_original)))]
# res_gwfo[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]
# 
# res_gwfo<-UnifyLoci(res_gwfo,variant_col = 'variant_ID',
#                    locus_col = 'locuscontext_id',group.by = 'chr')
# unique(res_gwfo$uni.locus_id)|>length()#174/195
# unique(res_gwfo[str_detect(credibleset,'cs95')]$uni.locus_id)|>length()#92/91
# 


#ADxQTL coloc
#run gwas_cs_extension.R
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#add the coloc locus overlap between studies
res_c[,chr:=seqid(variant_ID)]
res_c[,pos:=pos(variant_ID)]
res_c<-UnifyLoci(res_c[,-'uni.locus_id'],variant_col = 'variant_ID',
                 locus_col = 'locuscontext_id',group.by = c('chr'),rm_overlap = T)
unique(res_c$uni.locus_id)#146
table(res_c$context)
setdiff(res_c$context,contexts$context)
fwrite(res_c,fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

mtd[Method=='ColocBoost'&Modality=='AD_xQTL_colocalization',summary_file:=fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]
mtd[Method=='ColocBoost'&Modality=='AD_xQTL_colocalization']


#ADAD colocs
res_cad2<-fread(fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

mtd[Method=='ColocBoost'&Modality=='AD_meta_colocalization',summary_file:=fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]
mtd[summary_file=='']



#fsusie Coloc
res_cf<-fread(fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#add the coloc locus overlap between studies
res_cf[,chr:=seqid(variant_ID)]
res_cf[,pos:=pos(variant_ID)]
res_cf<-UnifyLoci(res_cf[,-'uni.locus_id'],variant_col = 'variant_ID',
                 locus_col = 'locuscontext_id',group.by = c('chr'),rm_overlap = T)
unique(res_cf$uni.locus_id)#24
table(res_cf$context)#
res_cf[context=='']#the extension
fwrite(res_cf,fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))
mtd[Method=='Coloc'][]
mtd[Method=='Coloc',summary_file:=fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz')]


#  4) the TWAS/MR####
table(mtd$Method)
res_twas<-fread(file.path('/data',mtd[Method=='twas']$Path))
setdiff(res_twas$context,res_scf$context)
res_twas[,gene_ID:=str_extract(molecular_id,'ENSG[0-9]+')]
res_twas[,event_ID:=context]
res_twas[str_detect(context,'_chr|_un'),context:=strsplit(context[1],'_chr|_gp_')[[1]][1],by='context']
setdiff(res_twas$context,res_scf$context)#ok

res_twas[,gwas_source:=paste0('AD_',gwas_study)]
res_twas<-res_twas[,-'gwas_study']

#TWAS signif
res_twas<-unique(res_twas[order(gene_ID,context,gwas_source,twas_pval)],by=c('gene_ID','context','gwas_source','method'))

#TWAS sig: (pvalue<2.5e-6 for >50% methods OR (pvalue<2.5e-6 AND the best method))
nmethods<-unique(res_twas$method)|>length()
res_twas[!is.na(twas_pval),TWAS_signif:=(twas_pval<2.5e-6&(is_selected_method))|(sum(twas_pval<2.5e-6)>=nmethods/2),
         by=c('gene_ID','block','context','gwas_source')]
res_twas[(TWAS_signif)]$gene_ID|>unique()|>length()#598 genes

fwrite(res_twas,fp(out,'res_AD_XWAS.csv.gz'))
res_twas<-fread(fp(out,'res_AD_XWAS.csv.gz'))

#add GWAS info: top gwas variant and combined locus for those variants

#get the xQTL CS with abs(GWAS z)>3
res_sc<-rbindlist(lapply(file.path('/data',mtd[Method=='single_context_finemapping']$Path),function(f)fread(f)),fill = T)
res_sc[,context:=str_extract(event_ID[1], '^.+?(?:_chr|_ENSG|_gp_)')%>%
         gsub('(_chr|_ENSG|_gp_)', '', .),by='event_ID']
setdiff(res_twas$context,res_sc$context)
#find all CS for those twas signif gene context
res_sc[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   ifelse(cs_coverage_0.5>0,paste0('cs50_',cs_coverage_0.5),'out_cs')))]
res_sc[,locuscontext_id:=paste(event_ID,credibleset,sep='_')]

res_twascs<-merge(res_twas,
                  unique(res_sc[order(locuscontext_id,-PIP,-abs(conditional_effect))][,.(context,gene_ID,
                                                                                         locuscontext_id,
                                                                                         variant_ID,PIP)],
                         by=c('locuscontext_id')),
                  by=c('context','gene_ID'),all.x = T,
                  allow.cartesian = TRUE)

#add gwas zscore for those CS
res_twascs<-rbindlist(lapply(unique(res_twascs$gwas_source),function(g){
  message(g)
  file=ps('/data/analysis_result/summary_stats_qced/AD_GWAS_RSS_QC_RAISS_imputed_concatenate_sumstat/',g,'_RSS_QC_RAISS_imputed.tsv.gz')
  if(file.exists(file)){
    gwas<-fread(file)
    merge(res_twascs[gwas_source==g],gwas[,.(variant_ID=paste0('chr',variant_id),
                                             gwas_zscore=z)],
          all.x = T,by='variant_ID')
  }else{
    res_twascs[gwas_source==g]
  }
  
}),fill = T)
#n twas gene with CS with gwas z >3
res_twascs[(TWAS_signif)]$gene_ID|>unique()|>length()#595
res_twascs[(TWAS_signif)][abs(gwas_zscore)>3]$gene_ID|>unique()|>length()#222/595
res_twascsf<-res_twascs[(TWAS_signif)][abs(gwas_zscore)>3]
summary(abs(res_twascsf$gwas_zscore))
fwrite(res_twascsf,fp(out,'res_AD_XWAS_topCS_GWASzscore3.csv.gz'))
res_twascsf<-fread(fp(out,'res_AD_XWAS_topCS_GWASzscore3.csv.gz'))
#QC: gene found in coloc but not in XWAS
setdiff(res_c$gene_ID,res_twascsf$gene_ID) #120/204
#comp to without filtering
setdiff(res_c$gene_ID,res_twas[(TWAS_signif)]$gene_ID) #92


#take top method only
res_twascsfu<-unique(res_twascsf[order(gene_ID,context,twas_pval,-abs(gwas_zscore))],
                    by=c('gene_ID','context','gwas_source'))

#add r0.8 twas union locus
library(pecotmr)
#first get long table variant1 variant2 r for abs(r)>0.5
#for all block
regions=unique(res_twascsfu$block)
variants_cors<-rbindlist(lapply(1:length(regions),function(i){
  message(i,'/',length(regions))
  b=regions[i]
  ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
                          region = data.frame(chr=seqid(b),
                                              start=start(b),
                                              end=end(b)))$combined_LD_matrix
  variants<-intersect(rownames(ldblock),str_remove(res_twascsfu$variant_ID,'^chr'))
  variants_dt<-data.table(variant_id=variants)
  variants_cors<-variants_dt[,{
    vars_cor<-rownames(ldblock)[abs(ldblock[,variant_id])>0.5]
    
    list(variant2=vars_cor,
         r=ldblock[vars_cor,variant_id])
  },by=.(variant_id)]
  variants_cors[,block:=b][]
}))
fwrite(variants_cors,fp(out,'twas_variants_cor0.5.csv.gz'))
variants_cors<-fread(fp(out,'twas_variants_cor0.5.csv.gz'))

variants_cors[,variant_ID:=paste0('chr',variant_id)]

variants_cors[,variant_ID2:=paste0('chr',variant2)]

#add those variants
res_twascsfu
res_twascsfu<-merge(res_twascsfu[,-'correlated_variants'],
                    unique(variants_cors[abs(r)>0.5][,.( correlated_variants=paste(variant_ID2,
                                                                                   collapse = '|')),
                                                     by='variant_ID']),
                   by='variant_ID')


res_twascsfu[,chr:=seqid(variant_ID)]
res_twascsfu[,gwas_pval:=getPval(gwas_zscore)]
res_twascsfu<-UnifyLoci(res_twascsfu[,-'uni.locus_id'],variant_col = 'variant_ID',
                         locus_col = 'locuscontext_id',
                         overlap_col = 'correlated_variants',
                         rm_overlap = TRUE,group.by = 'chr',
                        order.by = 'gwas_pval')
unique(res_twascsfu$uni.locus_id)#136/1k
res_twascsfu[uni.locus_id=='chr19_1']
fwrite(res_twascsfu,fp(out,'res_AD_XWAS_topCS_GWASzscore3_r0.5locus_union.csv.gz'))

#    - MR####
#export table: export_mr.R
res_mr<-fread('/adpelle1/export/FunGen_mr.exported.bed.gz')
res_mr[str_detect(context,'DLPFC_Bennett_pQTL'),context:='DLPFC_Bennett_pQTL']
res_mr[!is.na(Q_pval)][duplicated(paste(gene_ID,context,gwas_study))]

#harmonize gwas source /study
res_mr[,gwas_source:=paste0('AD_',gwas_study)]
table(res_mr$gwas_source)

# add to TWAS and call MR signif 
res_twas<-fread(fp(out,'res_AD_XWAS.csv.gz'))

res_mrtwas<-merge(unique(res_mr[order(gene_ID,context,gwas_source,-num_CS,I2)],by=c('gene_ID','context','gwas_source')),
                unique(res_twas[!is.na(twas_pval)][order(gene_ID,context,gwas_source,twas_pval)],by=c('gene_ID','context','gwas_source')),all=TRUE)
res_mrtwas[,MR_signif:=TWAS_signif&cpip>0.5&num_CS>=2&I2<0.5]

res_mrtwas[(TWAS_signif)]
res_mrtwas[(MR_signif)]
fwrite(res_mrtwas,fp(out,'res_AD_XWAS_MR.csv.gz'))
fwrite(res_mrtwas[(TWAS_signif)],fp(out,'res_AD_XWAS_MR_filtered_TWAS_sig.csv.gz'))

mtd[Method%in%c('twas','MR'),summary_file:=fp(out,'res_AD_XWAS_MR_filtered_TWAS_sig.csv.gz')]
fwrite(mtd,'metadata_analysis.csv')
unique(mtd[!is.na(summary_file)][,.(Method,summary_file)])
#confirm all have variant_ID, gene_ID, context, locuscontext_id
cols_to_check<-c('variant_ID', 'gene_ID', 'context','locuscontext_id')

unique(mtd[!is.na(summary_file)][,.(Method,summary_file)])[,{
  message(summary_file[1])
  missing<-setdiff(cols_to_check,colnames(fread(summary_file[1],nrows = 1)))
  if(length(missing)>0){
    warning('missing ',paste(missing,collapse = ' '))
    
  }else{
    message('OK')
  }
},by='summary_file']

#ok except TWAS/MR but normal
mtd[,variant_level_method:=!Method%in%c('MR','twas')]
fwrite(mtd,'metadata_analysis.csv')
mtd<-fread('metadata_analysis.csv')
mtd[summary_file=='']

#II) get unified 'AD potential' locus ####
#outputs:ADunilocus variant ADlocus_event ADmethod
#here we want to integrate GWAS locus found in i) single gwas finemapping, ii) ADxQTL coloc, iii) ADxAD colocs

res_gwf<-fread(mtd[Method=='AD_GWAS_finemapping']$summary_file[1])
unique(mtd$Modality)
res_cad<-fread(mtd[Modality=='AD_xQTL_colocalization']$summary_file)[str_detect(context,'^AD')]

res_cad2<-fread(mtd[Modality=='AD_meta_colocalization']$summary_file)
res_cad3<-fread(mtd[Method=='Coloc']$summary_file[1])

# #with TWAS locus?
# res_twascsfu<-fread(fp(out,'res_AD_XWAS_topCS_GWASzscore3_r0.5locus_union.csv.gz'))
# unique(res_twascsfu$uni.locus_id)
# #keep only top uni locus based on gwas z
# 
# res_twascsfuf<-unique(res_twascsfu[order(-abs(gwas_zscore))],by=c('uni.locus_id'))
# 
# res_ad<-rbindlist(list(res_gwf[,.(variant_ID,locuscontext_id,AD_method='single_gwas_finemapping',gwas_source=event_ID,credibleset,PIP)],
#                        res_c[,.(variant_ID,locuscontext_id,AD_method='ADxQTL_coloc',gwas_source=gwas_source,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
#                        res_cad2[,.(variant_ID,locuscontext_id,AD_method='ADxAD_coloc',gwas_source=event_ID,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
#                        res_twascsfuf[,.(variant_ID,locuscontext_id,twas_unilocus=paste0('TWAS_',uni.locus_id),gwas_source,
#                                         AD_method='TWAS',twas_z)]),fill = T)

res_ad<-rbindlist(list(res_gwf[,.(variant_ID,locuscontext_id,AD_method='single_gwas_finemapping',gwas_source=event_ID,credibleset,PIP)],
                       res_cad[,.(variant_ID,locuscontext_id,AD_method='ADxQTL_coloc',gwas_source=gwas_source,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad2[,.(variant_ID,locuscontext_id,AD_method='ADxAD_coloc',gwas_source=event_ID,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad3[,.(variant_ID,locuscontext_id,AD_method='ADFsusie_coloc',gwas_source=gwas_source,SNP_PPH4,L_PP.H4.abf,xQTL_L,AD_L)]),fill = T)
res_ad[,chr:=seqid(variant_ID)]


res_ad<-UnifyLoci(res_ad,variant_col = 'variant_ID',locus_col = 'locuscontext_id',
                  group.by = 'chr')
setnames(res_ad,'uni.locus_id','ADlocus')

# #twas spe
# res_ad[,n.meth.locus:=length(unique(AD_method)),by='ADlocus']
# unique(res_ad[n.meth.locus==1],by='ADlocus')$AD_method|>table()
# 
#             # ADxAD_coloc            ADxQTL_coloc single_gwas_finemapping                    TWAS 
#             #           8                      63                      70                      40 
# 
# #for TWAS only locus in chr19, keep only the top hit
# res_ad<-res_ad[!(n.meth.locus==1&AD_method=='TWAS'&chr=='chr19'&twas_unilocus!='TWAS_chr19_1')]
# unique(res_ad$ADlocus)#289
# 
# unique(res_ad[n.meth.locus==1],by='ADlocus')$AD_method|>table()
#           # ADxAD_coloc            ADxQTL_coloc single_gwas_finemapping                    TWAS 
#           #             8                      50                      53                      40 


#add gwas zscore
res_ad<-rbindlist(lapply(unique(res_ad$gwas_source),function(g){
  message(g)
  file=ps('/data/analysis_result/summary_stats_qced/AD_GWAS_RSS_QC_RAISS_imputed_concatenate_sumstat/',g,'_RSS_QC_RAISS_imputed.tsv.gz')
  if(file.exists(file)){
    gwas<-fread(file)
    merge(res_ad[gwas_source==g],gwas[,.(variant_ID=paste0('chr',variant_id),
                                             gwas_zscore=z)],
          all.x = T,by='variant_ID')
  }else{
    res_ad[gwas_source==g]
  }
  
}),fill = T)

#ncs95 or ADxAD coloc
res_ad[str_detect(credibleset,'cs95')]$ADlocus|>unique()|>length()#88
res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc']$ADlocus|>unique()|>length()#116
res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc'|AD_method=='ADxAD_coloc'|AD_method=='ADxQTL_coloc']$ADlocus|>unique()|>length()#191
#res_ad[str_detect(credibleset,'cs95')|AD_method=='ADxAD_coloc'|AD_method=='ADxAD_coloc'|AD_method=='ADxQTL_coloc'|AD_method=='TWAS']$ADlocus|>unique()|>length()#191
res_ad$ADlocus|>unique()|>length()#249

#add more meaningful ADLocus name: chr-[start]-[end]
res_ad[,is.cs95:=any(str_detect(credibleset,'cs95'),na.rm = T),by='ADlocus']
res_ad[,ADlocusID:=paste(seqid(ADlocus),min(pos(variant_ID)),max(pos(variant_ID)),sep = '_'),by='ADlocus']

res_ad[,ADlocus.n.variant.union:=length(unique(variant_ID)),by='ADlocus']
setcolorder(res_ad,c('ADlocusID','ADlocus.n.variant.union','is.cs95'))

table(unique(res_ad,by=c('ADlocus','variant_ID'))$ADlocusID)

#save
unique(res_ad[(is.cs95)]$ADlocusID)
fwrite(res_ad,fp(out,'AD_loci_unified2.csv.gz'))
res_ad<-fread(fp(out,'AD_loci_unified2.csv.gz'))
res_ad[ADlocus=='chr1_1']
fread(mtd[Modality=='AD_xQTL_colocalization']$summary_file)[locuscontext_id=='ENSG00000288636:cos1:y1_y5_AD_Bellenguez_EADB_2022_ENSG00000288636']


unique(res_ad[!is.na(coloc_coef)][abs(coloc_coef)<1e-5&sign(gwas_zscore)!=sign(coloc_coef)][,paste(variant_ID,gwas_source)])|>length()
unique(res_ad[!is.na(coloc_coef)][abs(coloc_coef)<1e-5][,paste(variant_ID,gwas_source)])|>length()#267/4602

unique(res_ad[!is.na(coloc_coef)][,paste(gwas_source,AD_method)])
table(res_ad[!is.na(coloc_coef)][sign(gwas_zscore)!=sign(coloc_coef)][,paste(gwas_source,AD_method)])

unique(res_ad[abs(coloc_coef)>1e-5&sign(gwas_zscore)!=sign(coloc_coef)][,paste(variant_ID,gwas_source)])|>length()#2123
unique(res_ad[abs(coloc_coef)>1e-5][,paste(variant_ID,gwas_source)])|>length()#2123/10837

#top variants per cos id
res_ad[,top_vcp:=vcp==max(vcp,na.rm = T),by=c('ADlocus','gwas_source')]
unique(res_ad[(top_vcp)][!is.na(coloc_coef)][abs(coloc_coef)<1e-5&sign(gwas_zscore)!=sign(coloc_coef)][,paste(variant_ID,gwas_source)])|>length()
unique(res_ad[(top_vcp)][!is.na(coloc_coef)][abs(coloc_coef)<1e-5][,paste(variant_ID,gwas_source)])|>length()#267/4602

unique(res_ad[(top_vcp)][abs(coloc_coef)>1e-5&sign(gwas_zscore)!=sign(coloc_coef)][,paste(variant_ID,gwas_source)])|>length()#2123
unique(res_ad[(top_vcp)][abs(coloc_coef)>1e-5][,paste(variant_ID,gwas_source)])|>length()#24/207
# #because summary qced mismatch exported?
# res_sc<-fread('analyses_summary/res_all_single_gwas_finemapping_cs50orgreater_extended_unified.csv.gz')
# res_scf<-res_sc[!is.na(PIP)]
# res<-merge(res_scf[,.(variant_ID,gwas_source,z)],unique(res_ad[,.(variant_ID,gwas_source,gwas_zscore)]),by=c('variant_ID','gwas_source'))

# res[sign(gwas_zscore)!=sign(z)]

#one variant per line
#with  chr, pos, effect_allele, GWAS methods,
res_ad[,variant_inclusion_probability:=ifelse(is.na(vcp),PIP,vcp)]
res_ad[is.na(variant_inclusion_probability),variant_inclusion_probability:=SNP_PPH4]

res_ad[,max_variant_inclusion_probability:=max(variant_inclusion_probability,na.rm = T),by='variant_ID']

res_ad[,max_variant_inclusion_probability_method:=AD_method[which.max(variant_inclusion_probability)],by='variant_ID']
res_ad[,n.variants:=length(unique(variant_ID)),by='ADlocus']

res_ad[,chr:=seqid(variant_ID,only_num =T )]
res_ad[,pos:=pos(variant_ID)]
res_ad[,max_zscore:=gwas_zscore[which.max(abs(gwas_zscore))],by='variant_ID']
res_ad[,min_pval:=getPval(max_zscore)]
res_ad[,susie_coverage:=sort(unique(str_extract(credibleset,'cs[0-9]+')),decreasing = T)[1],by='ADlocusID']

res_ad[order(variant_ID,-variant_inclusion_probability),GWAS_methods:=paste(unique(AD_method),collapse =  '|'),by='variant_ID']

res_ad[order(variant_ID,-variant_inclusion_probability),gwas_sources:=paste(unique(gwas_source),collapse =  '|'),by='variant_ID']
#res_ad[,gwas_source_effect:=NULL]
res_ad[!is.na(gwas_zscore),gwas_source_effect:=paste0(str_extract(gwas_source,'[BWJK]'),ifelse(gwas_zscore>0,'+','-'))]
res_ad[is.na(gwas_zscore),gwas_source_effect:=str_extract(gwas_source,'[BWJK]')]

#res_ad[!is.na(coloc_coef)&!is.na(gwas_zscore),gwas_source_effect:=paste0(gwas_source_effect,ifelse(coloc_coef>0,'+','-'))]
# res_ad[!is.na(coloc_coef)&!is.na(gwas_zscore),gwas_source_effect.locus:=gwas_source_effect[order(-str_length(gwas_source_effect))][1],by=.(variant_ID,gwas_source)]

res_ad[order(variant_ID,-variant_inclusion_probability),gwas_sources_effects:=paste(gwas_source_effect[!is.na(gwas_source_effect)&!duplicated(str_extract(gwas_source_effect,'[BWJK]'))],collapse =  '|'),by='variant_ID']
unique(res_ad[,.(ADlocusID,gwas_sources_effects)])[1:100]


res_ad[!is.na(coloc_coef)&!is.na(gwas_zscore)][sign(gwas_zscore)!=sign(coloc_coef)]
#remove variants with min_pval>1e-5 or in cs70 or cs50 only
res_adf<-res_ad[min_pval<1e-5]
unique(res_adf$ADlocusID) #203/249
#only cs95 or coloc
res_adf<-res_adf[ADlocusID%in%ADlocusID[susie_coverage=='cs95'|str_detect(GWAS_methods,'coloc')]]
length(unique(res_adf$ADlocusID)) #197/249
unique(res_adf[min_pval<5e-8]$ADlocusID)|>length() #103
unique(res_adf[susie_coverage=='cs95']$ADlocusID)|>length()  #88
res_adf[is.na(is.cs95),is.cs95:=FALSE]


#how the single method loci of those 166 are far away from the closest? 
#TODO

res_adfv<-unique(res_adf[order(variant_ID,-abs(gwas_zscore),-variant_inclusion_probability)][,.(ADlocusID,
                                                                                             n.variants,
                                                                                          chr,
                                                                                          pos,
                                                                                          variant_ID,
                                                                                          GWAS_methods,
                                                                                          max_variant_inclusion_probability,
                                                                                          max_variant_inclusion_probability_method,
                                                                                          gwas_sources,
                                                                                          gwas_sources_effects,
                                                                                          is.cs95,
                                                                                          max_zscore,
                                                                                          min_pval,
                                                                                          susie_coverage,
                                                                                          ADlocus)],
                by='variant_ID')

res_adfv[,max_variant_inclusion_probability_rank:=rank(-max_variant_inclusion_probability),by='ADlocus']
res_adfv[,n.variants_0.1:=sum(max_variant_inclusion_probability>0.1),by='ADlocus']

res_adfv[order(chr,ADlocus,-max_variant_inclusion_probability)]
res_adfv[is.na(susie_coverage)]
table(res_adfv$ADlocusID)
#unique(res_adfv[,.(ADlocus,GWAS_methods)][order(-str_length(GWAS_methods))],by='ADlocus')[101:194]

#add locus index
adloci<-unique(res_adfv[order(chr,ADlocus,-max_variant_inclusion_probability)],by='ADlocus')
adloci[,locus_index:=1:.N]
res_adfv<-merge(res_adfv[,-'locus_index'],adloci[,.(ADlocus,locus_index)])
setcolorder(res_adfv,'locus_index')
unique(res_adfv[str_detect(GWAS_methods,'Fsusie')],by='ADlocusID')|>nrow()#18
unique(res_adfv[GWAS_methods=='ADFsusie_coloc'],by='ADlocusID')|>nrow()#2

#add start and end


#add rsid
#FIXME ref contain only chr start end to identify the variant, lacking REF/ALT
# variants<-fread('/data/resource/references/00-All.variants.gz')
# setnames(variants,c('chr','start','end','rsid'))
# res_adfv<-merge(res_adfv,variants[chr%in%as.character(1:22)][,.(chr=as.numeric(chr),pos=end,rsid)],
#                 by=c('chr','pos'),all.x=TRUE)
# res_adfv[is.na(rsid)]#528/6100
# 
# res_adfv<-rbind(res_adfv[!is.na(rsid)],merge(res_adfv[is.na(rsid)][,-'rsid'],variants[chr%in%as.character(1:22)][,.(chr=as.numeric(chr),pos=start,
#                                                                              rsid)],by=c('chr','pos'),all.x=TRUE))
# 
# res_adfv[is.na(rsid)]#528/6100
# res_adfv<-rbind(res_adfv[!is.na(rsid)],
#                 merge(res_adfv[is.na(rsid)][,-'rsid'],variants[chr%in%as.character(1:22)][,.(chr=as.numeric(chr),pos=start+1,
#                                                                                              rsid)],
#                       by=c('chr','pos'),all.x=TRUE))
# 
# res_adfv<-rbind(res_adfv[!is.na(rsid)],
#                 merge(res_adfv[is.na(rsid)][,-'rsid'],variants[chr%in%as.character(1:22)][,.(chr=as.numeric(chr),pos=end+1,
#                                                                                                           rsid)],
#                        by=c('chr','pos'),all.x=TRUE))
# res_adfv[is.na(rsid)]#10/6100
# 
# 
# variants[end-start>1]

#add v2f score
files<-list.files('/data/analysis_result/cv2f/score/',pattern = '.cv2f$',full.names = TRUE)
v2f<-rbindlist(lapply(files,function(f)fread(f)))
v2f[!str_detect(SNP,'rs')]
res_adfv<-merge(res_adfv,v2f[,.(chr=CHR,pos=BP,cV2F,SNP)],
                      by=c('chr','pos'),all.x=TRUE)


res_adfv[is.na(cV2F)]
res_adfv[!is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]
res_adfv[is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]
res_adfv[is.na(cV2F)][str_count(alt(variant_ID))>1|str_count(ref(variant_ID))>1]$max_variant_inclusion_probability|>summary()


res_adfv[!is.na(cV2F),cV2F_rank:=rank(-cV2F),by='ADlocus']

fwrite(res_adfv[order(locus_index)],fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))
res_adfv<-fread(fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))


#export
fwrite(res_adfv[order(locus_index)],fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))




#III) Get the merged long variant-level table  AD xQTL overlap: each row a variant-ADlocus-Method-context-gene####
#add those ADlocus annot to all variant level summ table
mtd<-fread('metadata_analysis.csv')
res_adfv<-fread('analyses_summary/AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz')
update_summary_ad=FALSE

mtd[file.exists(summary_file)]

#first check if need flipping to harmonize variant_ID
mismatches<-unique(mtd,by='summary_file')[file.exists(summary_file)&(variant_level_method),{
  message(Method[1])

    res<-fread(summary_file[1])
    res[,chr:=seqid(variant_ID,only_num = T)]
    res[,pos:=pos(variant_ID)]
    res[,variant_id:=variant_ID]
    
    res<-merge(res[,-c('ADlocus','locus_index','ADlocusID','variant_ID')],
               res_adfv[,.(variant_ID,chr,pos,ADlocus,locus_index,ADlocusID)],all.x=T,by=c('chr','pos'))
    print(table(res[!is.na(variant_ID)][variant_id!=variant_ID]$context))
    res[!is.na(variant_ID)][variant_id!=variant_ID][,.(variant_id,variant_ID,context)]
},by=c('Method','Modality')]
unique(mismatches$variant_ID)#213/5318 variants mismatch
unique(mismatches[,.(variant_id,variant_ID)])

effect_cols<-c('conditional_effect','coef','effect_direction','top_effect')

mismatches[,to_flip:=ref(variant_id)==alt(variant_ID)|alt(variant_id)==ref(variant_ID)]
mismatches[(to_flip)]
table(mismatches[(to_flip)]$context)
table(mismatches[(to_flip)]$Method)
table(mismatches[(to_flip)][,paste(Method,context)])
fwrite(mismatches,fp(out,'variant_mismatchs_with_AD2.csv.gz'))

mismatches<-fread(fp(out,'variant_mismatchs_with_AD2.csv.gz'))


#Annotate each variant level Methods table  with the AD Loci
mtd[file.exists(summary_file)&(variant_level_method),summary_file_ad:={
    message(Method[1])
    summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file[1],compression = T),'_overlapADloci.csv.gz')
    
    if(!file.exists(summary_file_ad)|update_summary_ad){
      res<-fread(summary_file[1])
      res[,chr:=seqid(variant_ID,only_num = T)]
      res[,pos:=pos(variant_ID)]
      res[,variant_id:=variant_ID]
      
      if('locuscontext_id'%in%colnames(res)){
        #n.variants
        res[,n.variant.locus:=.N,by='locuscontext_id']
      }
      res<-merge(res[,-c('ADlocus','locus_index','ADlocusID','variant_ID')],
                 res_adfv[,.(variant_ID,chr,pos,ADlocus,locus_index,ADlocusID)],
                 all.x=T,by=c('chr','pos'))
      
      res[,is.in.ad.locus:=!is.na(ADlocus)]
      
      res[is.na(variant_ID),variant_ID:=variant_id]
      
      if(nrow(mismatches)>0){
        #flip effect if needed
        effect_colsf<-intersect(effect_cols,colnames(res))
        if(length(effect_colsf)>0){
          variants_to_flip<-intersect(res$variant_id,mismatches$variant_id)
          
          if(length(variants_to_flip)>0){
            res<-rbind(res[variant_id%in%variants_to_flip,(effect_colsf):=lapply(.SD,function(x)-x),.SDcols = effect_colsf],
                       res[!variant_id%in%variants_to_flip])
          }
          
        }
      }
      
      
      
      if('locuscontext_id'%in%colnames(res)){
        
        res[,n.variant.overlap.ad:=sum(is.in.ad.locus),by=.(locuscontext_id)]
        res[,overlap.ad:=any(is.in.ad.locus),by=.(locuscontext_id)]
        resf<-res[(overlap.ad)]
        
      }else{
        resf<-res[(is.in.ad.locus)]
        
      }
      
      fwrite(res,summary_file[1])
      
      
      fwrite(resf,summary_file_ad)
    }
    
    summary_file_ad
  },by=c('summary_file')]


fwrite(mtd,'metadata_analysis.csv')
mtd<-fread('metadata_analysis.csv')


#get the long AD overlap table
res_adx<-rbindlist(lapply(unique(mtd[file.exists(summary_file_ad)&(variant_level_method)]$summary_file_ad),function(f){
  message(f)
  m=mtd[summary_file_ad==f]$Method[1]
  
  if(m=='ColocBoost'){
    m=mtd[summary_file_ad==f]$Modality[1]
    summary_file_ad=mtd[Modality==m]$summary_file_ad[1]
    
  }else{
    summary_file_ad=mtd[Method==m]$summary_file_ad[1]
    
  }
  message(summary_file_ad)
  res<-fread(summary_file_ad)[,Method:=m]
  message('found associated with ',length(unique(res$ADlocus)),' AD locus')
  return(res)
}),fill = T)
res_adx
res_adx[!str_detect(gene_ID,'^ENSG')]$context|>table() #OK



#populate all locuscontext_id lvl with the AD locus overlapping it
res_adx[ADlocus=='',ADlocus:=NA]
res_adx[,ADlocus:=ifelse(is.na(ADlocus),ADlocus[order(-n.variant.overlap.ad,is.na(ADlocus))][1],ADlocus),
    by='locuscontext_id']
res_adx[ADlocusID=='',ADlocusID:=NA]
res_adx[,ADlocusID:=ifelse(is.na(ADlocusID),ADlocusID[order(-n.variant.overlap.ad,is.na(ADlocusID))][1],ADlocusID),
    by='locuscontext_id']
res_adx[,locus_index:=ifelse(is.na(locus_index),locus_index[order(-n.variant.overlap.ad,is.na(locus_index))][1],locus_index),
    by='locuscontext_id']

#rm variants from extension in the AD methods
res_adx<-res_adx[!(is.na(vcp)&str_detect(Method,'coloc'))]
res_adx<-res_adx[!(is.na(PIP)&str_detect(Method,'finemapping'))]

#add GWAS infos : 
dupcols<-setdiff(intersect(colnames(res_adfv),colnames(res_adx)),'variant_ID')
res_adx<-merge(res_adx[,.SD,.SDcols = !dupcols],res_adfv,
               by=c('variant_ID'),all = T)
res_adx[locus_index==1]

# #add gwaszscore for gwas method
# res_adx<-merge(res_adx,unique(res_ad[,.(variant_ID,context=gwas_source,gwas_zscore)]),all.x = T,by=c('variant_ID','context'))
res_adx[(is.in.ad.locus)]|>nrow()

#bind to TWAS/MR
#for twas, get ADlocus-gene connection to annotate 
summary_file=mtd[Method=='twas']$summary_file[1]
summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file,compression = T),'_overlapADloci.csv.gz')
update_twasad<-TRUE
if(!file.exists(summary_file_ad)|update_twasad){
  res_adgenes<-unique(res_adx[,.(ADlocus,gene_ID,ADlocusID,locus_index)])
  
  
  res_tw<-fread(summary_file)
  res_twad<-merge(res_tw[,-c('ADlocus','ADlocusID','locus_index')],res_adgenes,by='gene_ID',all.x = T,allow.cartesian = T)
  res_twadf<-res_twad[!is.na(ADlocus)][(TWAS_signif)]
  unique(res_twadf$gene_ID)
  fwrite(res_twadf,summary_file_ad)
  mtd[Method%in%c('twas','MR'),summary_file_ad:=..summary_file_ad]
  fwrite(mtd,'metadata_analysis.csv')
}


res_twadf<-fread(summary_file_ad)
res_twadf[(MR_signif)]
res_twadf[is.na(MR_signif),MR_signif:=FALSE]

#bind it
res_adx<-rbind(res_adx,res_twadf[,Method:='TWAS/MR'],fill=T)

#populate TWAS info at variant level
res_adx[,TWAS_signif:=any(TWAS_signif),by=.(ADlocus,gene_ID,context,gwas_source)]
res_adx[,MR_signif:=any(MR_signif),by=.(ADlocus,gene_ID,context,gwas_source)]


table(unique(res_adx,by=c('locuscontext_id','ADlocus','gene_ID','context'))$Method)
# AD_GWAS_finemapping     AD_meta_colocalization     AD_xQTL_colocalization           APOE interaction 
# 341                        249                       2458                         39 
# Coloc         fSuSiE_finemapping           msex interaction  multi_context_finemapping 
# 93                        137                          2                        747 
# multi_gene_finemapping                         QR single_context_finemapping          trans_finemapping 
# 527                        848                       2333                          2 
# TWAS/MR 
# 661 


#add context broad and short
#Change FOR sQTL: if 'UP' put as u-sQTL (Unproductive sQTL), if 'PR' as p-sQTL (productive)
res_adx[str_detect(event_ID,':UP:')&!str_detect(context,'u_sQTL'),context:=str_replace(context,'sQTL','u_sQTL')]
res_adx[str_detect(event_ID,':PR:')&!str_detect(context,'p_sQTL'),context:=str_replace(context,'sQTL','p_sQTL')]
res_adx[str_detect(event_ID,':UP:')]$context

contexts<-fread('contexts_metadata.csv')
setdiff(res_adx$context,contexts$context) 
setdiff(contexts$context,res_adx$context)

res_adx<-merge(res_adx[,-c('context_broad','context_short')],
               unique(contexts[,.(context,context_broad,context_short)]),all.x = T,by='context')
table(res_adx[,.(Method,context)])

table(res_adx[,.(Method,context_short)])


#add lacking contexts
res_adx[is.na(context)]$locuscontext_id
res_adx[is.na(context),context:=str_extract(locuscontext_id,'metabolome')]


#GENE info tss /tes
update_geneinfo=FALSE
if(update_geneinfo){
  trans<-fread('/data/resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list')
  
  gtf<-fread('/data/resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf',
             select = c(1,3,4,5,7,9),col.names = c('chr','biotype','start','end','sens','info'))
  gtf[,gene_id:=str_extract(info,'ENSG[0-9]+')]
  trans<-merge(trans,gtf[biotype=='gene'][,.(gene_id,tss=ifelse(sens=='+',start,end),tes=ifelse(sens=='+',end,start))])
  setnames(trans,'gene_id','gene_ID')
  trans<-unique(trans,by='gene_ID')
  fwrite(trans,fp(out,'genes_infos.csv.gz'))
}else{
  trans<-fread(fp(out,'genes_infos.csv.gz'))
}

res_adx<-merge(res_adx[,-c('tss','tes','gene_name')],trans[,-c('start','end','#chr')],all.x = T,by='gene_ID')

res_adx[,chr:=seqid(variant_ID[1],only_num = T),by='variant_ID']
res_adx[,pos:=pos(variant_ID[1]),by='variant_ID']

res_adx[,effect_allele:=alt(variant_ID[1]),by='variant_ID']
res_adx[,distance_from_tss:=pos-tss]
res_adx[,distance_from_tes:=pos-tes]

#associate to a gene each epiQTL: populate gene name for epiQTL if geneQTL overlap the locus
xqtl_methods<-c('fSuSiE_finemapping','single_context_finemapping',
                'multi_context_finemapping','AD_xQTL_colocalization','TWAS/MR')

res_adx<-rbind(res_adx[Method!='fSuSiE_finemapping'|is.na(Method)],
               res_adx[Method=='fSuSiE_finemapping']|>split(by=c('locus_index','context'))|>lapply(function(epi){
                 epivars<-epi$variant_ID
                 genes_overlapping_epi<-res_adx[Method%in%setdiff(xqtl_methods,'fSuSiE_finemapping')][variant_ID%in%epivars]$gene_name|>unique()
                 
                 epi_genes<-rbindlist(lapply(genes_overlapping_epi,function(g)epi[,gene_name:=g]))
                 return(epi_genes)
               })|>rbindlist())
res_adx[,gene_ID:=gene_ID[!is.na(gene_ID)][1],by='gene_name']

res_adx[locus_index==1]


#xQTL level variant prioritization
#MAX_VIP for xQTL like for GWAS 
#i.e. max_vip, vip_rank, top_method, top_context pER GENE

res_adx[,variant_inclusion_probability_xqtl:=max(ifelse(is.na(vcp),PIP,vcp)[Method%in%xqtl_methods&!str_detect(context,'^AD')],na.rm = T),by=c('variant_ID','gene_ID')]
res_adx[is.infinite(variant_inclusion_probability_xqtl),variant_inclusion_probability_xqtl:=NA]
res_adx[is.na(variant_inclusion_probability_xqtl)] #those not asso to xQTL


#add variant rank
res_adxv<-unique(res_adx,by=c('variant_ID','gene_ID'))
res_adxv[,variant_rank_xqtl:=rank(-variant_inclusion_probability_xqtl),
         by=c('ADlocus','gene_ID')]
res_adx<-merge(res_adx[,-'variant_rank_xqtl'],res_adxv[,.(variant_ID,ADlocus,gene_ID,variant_rank_xqtl)],by=c('ADlocus','gene_ID','variant_ID'),all.x=T)
res_adx[,.(locuscontext_id,variant_inclusion_probability_xqtl,variant_rank_xqtl)]
#add top context associated
res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD'),method_xqtl:=Method[which.max(ifelse(is.na(vcp),PIP,vcp))],by=c('variant_ID','gene_ID')]
res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD'),context_xqtl:=context[which.max(ifelse(is.na(vcp),PIP,vcp))],by=c('variant_ID','gene_ID')]
res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD'),n.variant_xqtl:=n.variant.locus[which.max(ifelse(is.na(vcp),PIP,vcp))],by=c('variant_ID','gene_ID')]
res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD'),effect_xqtl:=ifelse(str_detect(method_xqtl,'fSuSiE'),top_effect,ifelse(str_detect(method_xqtl,'finemapping'),conditional_effect,coef)),by=c('variant_ID','gene_ID')]
res_adx[Method%in%xqtl_methods&!str_detect(context,'^AD'),coverage_xqtl:=str_extract(credibleset,'cs[0-9]+')[which.max(ifelse(is.na(vcp),PIP,vcp))],by=c('variant_ID','gene_ID')]
colnames(res_adx)

res_adx[locus_index==1]
#add broader context and qtl type
res_adx[,context_broad2:=ifelse(context_broad%in%c('bulk_monocyte_eQTL','bulk_macrophage_eQTL','bulk_microglia_eQTL'),'Immune (bulk)',
                                ifelse(str_detect(context_broad,'bulk_brain'),'Brain',
                                       str_remove(context_broad,'_eQTL|_snATAC|_brain')))]
res_adx[,qtl_type:=str_extract(context_broad,'(u_|p_)?[a-z]+QTL')]

res_adx[str_detect(context_broad,'snATAC'),qtl_type:='caQTL']

fwrite(res_adx,fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))

#keep only columns of interest
# #filter cols: keep only AD locus, gene, method, context, variant_id and gwas infos,  and key method specific infos
colnames(res_adx)|>cat(sep='\n')

cols<-fread('all_analysis_long_selected_columns.csv')$column_name
setdiff(cols,colnames(res_adx))
setdiff(colnames(res_adx),cols)

colsf<-intersect(cols,colnames(res_adx))
res_adxf<-res_adx[,..colsf]
res_adxf<-res_adxf[order(locus_index)]
res_adxf

fwrite(res_adxf,fp(out,'res_allanalysis_ADloci_overlap_selectedcols.csv.gz'))
res_adxf<-fread(fp(out,'res_allanalysis_ADloci_overlap_selectedcols.csv.gz'))

# locus level summary
# #outputs: ADlocus ADlocus_event ADmethod ADvariants gene context locus_event method n_variants others_info (n_overlap, pip, TWAS_Z..)
# #res_adxf[,leading_gwas_variant:=variant_ID[order(-max_variant_inclusion_probability)][1],by='ADlocus']
# res_adxl<-unique(res_adxf[order(-max_variant_inclusion_probability,-PIP,-vcp,-abs(max_zscore))],by=c('ADlocus',
#                                                                    'locuscontext_id',
#                                                                    'gene_ID',
#                                                                    'context','Method'))
# 
# 
# 
# 
# 
# res_adxl<-res_adxl[order(locus_index)]
# res_adxl
# fwrite(res_adxl,fp(out,'res_allanalysis_ADloci_overlap_locus_summary.csv'))

# res_adxl$locus_index|>unique()|>length()
# 
# #genome wide signif only
# res_adxlf1<-res_adxl[leading_gwas_variant_pval<5e-8]
# fwrite(res_adxlf1,fp(out,'res_allanalysis_ADloci_overlap_locus_summary_GWASpval_genomewide_signif.csv'))
# 
# #cs95 only
# res_adxlf2<-res_adxl[susie_coverage=='cs95']
# fwrite(res_adxlf2,fp(out,'res_allanalysis_ADloci_overlap_locus_summary_GWASpval_cs95_only.csv'))



#   AD loci sanity check####
#in coloc only, how many from one context only?
# gwas_methods<-c('AD_GWAS_finemapping','AD_meta_colocalization',
#                 'AD_xQTL_colocalization','Coloc')
# res_adxa<-res_adx[Method%in%gwas_methods]
# res_adxa[,ADxQTLcoloc_only:=all(Method=='AD_xQTL_colocalization')|all(Method=='Coloc'),by='ADlocus']
# res_adxa<-res_adxa[(ADxQTLcoloc_only)][Method=='AD_xQTL_colocalization'][!str_detect(context,'^AD')]
# table(res_adxa$context)
# 
# res_adxa[,onecontext_only:=length(unique(context[!str_detect(context,'^AD')]))==1,by='ADlocus']
# ggplot(unique(res_adxa,by='ADlocus'))+
#   geom_bar(aes(x=context,fill=onecontext_only))+
#   scale_x_discrete(guide=guide_axis(angle=60))+ggtitle('Coloc specific AD locus')+labs(fill='found only in this context')
# 
# res_adx[locus_index==74]$context|>unique()
# res_adx[locus_index%in%c(73,74,75)][,.(locus_index,chr,pos)][order(locus_index,pos)]|>unique(by='locus_index')
# 
# res_adx[context_broad=='bulk_brain_sQTL'][str_detect(Method,'coloc')]$event_ID
# 
# #stats loci
# #n loci (+\- xQTL annotated)
# #n genome wide signif (+\- xQTL annotated)
# # cs95 (+\- xQTL annotated)
# res_adx[,xqtl_annotated_locus:=any(!is.na(context)&context!=''&!str_detect(context,'^AD')),by='ADlocus']
# stats<-rbindlist(list(unique(res_adx[order(-xqtl_annotated_locus)],by='ADlocus')[,stringence:='p<1*10-5'],
#                       unique(res_adx[min_pval<5e-8][order(-xqtl_annotated_locus)][,stringence:='p<5*10-8'],by='ADlocus'),
#                       unique(res_adx[susie_coverage=='cs95'][order(-xqtl_annotated_locus)][,stringence:='singlegwas_susie_rss_cs95'],
#                              by='ADlocus')))
# 
# 
# ggplot(stats)+geom_bar(aes(x=stringence,fill=xqtl_annotated_locus))+labs(y='Number of AD Loci')
# 
# #distrib distance betweem loci
# #max pos locus n minus min pos locus n+1
# dists<-unique(res_adx[order(locus_index,pos)],by=c('ADlocus','pos'))[,.(locus_index,
#                                                                         distance_to_next=sapply(locus_index,function(i){
#   max_pos<-max(pos[locus_index==i])
#   min_pos_next<-min(pos[locus_index==i+1])
#   return(min_pos_next-max_pos)
#   
# })),by='chr']|>unique()
# ggplot(dists[distance_to_next!=-Inf])+
#   geom_histogram(aes(x=distance_to_next))+scale_x_log10()
# 
# dists[distance_to_next<1e6]$locus_index|>unique()|>length()
# dists$locus_index|>unique()|>length()




#IV) WIDE TABLE CREATION  ####
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))


res_adxub<-WideTable(res_adx,group.by = 'context_short')

fwrite(res_adxub,fp(out,'res_AD_variants_xQTL.csv.gz'))
res_adxub<-fread(fp(out,'res_AD_variants_xQTL.csv.gz'))


#FILTER: keep only variants with GWAS PIP/VCP > 0.1, for maximum of 5
# if non of the variant has GWAS PIP/VCP > 0.1 we just show top one based on GWAS PIP/VCP and on xQTLPIP/VCP
res_adxub[,top_variants:=((max_variant_inclusion_probability>=0.1)&(max_variant_inclusion_probability_rank<=5|cV2F_rank<=5))|max_variant_inclusion_probability_rank==1|variant_rank_xqtl==1|cV2F_rank==1,by='ADlocus']
res_adxubf<-res_adxub[(top_variants)]
nrow(res_adxubf)#2517
unique(res_adxubf$gene_name)#444
unique(res_adxubf[min_pval<5e-8]$gene_name)#231


#add supplemental cols
res_adxubf[,gwas_significance:=ifelse(min_pval<5e-8,'genome wide',
                                        ifelse(min_pval<1e-6,
                                               'suggestive',
                                               'ns'))]
res_adxubf[,mlog10pval:=-log10(min_pval)]


#variant inclusion top confidence level
res_adxubf[,top_confidence:=str_extract(xQTL_effects,'C[0-9]')]
res_adxubf$top_confidence|>unique()
res_adxubf[(MR_signif)][,.(locus_index,gene_name)][!is.na(locus_index)]|>unique()
#Main Excel Sheet creation ####
#get the columns metadata ready

PrepColsMtd<-function(cols,colsmtd,res,qtl_orders=c('')){
  #replicate columns for each $variable adding the varnames
  variable_gp<-str_extract(cols$grandparent_column,'\\$[A-Za-z0-9_.]+')|>na.omit()|>str_remove('\\$')|>as.vector()
  
  variable_p<-str_extract(cols$parent_column,'\\$[A-Za-z0-9_.]+')|>na.omit()|>str_remove('\\$')|>as.vector()

  variable_gp_vec<-colsmtd[wildcard==variable_gp]$r_name
  
  variable_p_vec<-colsmtd[wildcard==variable_p]$r_name
  
  cols_extended<-rbindlist(lapply(variable_gp_vec,function(v1){
    color=colsmtd[wildcard==variable_gp][r_name==v1]$color
    grandparentname<-colsmtd[wildcard==variable_gp][r_name==v1]$excel_name
    
    cols_ext<-rbindlist(lapply(variable_p_vec,function(v2){
      name=paste(v1,v2,sep = '.')
      
   return(cols[expanded==1][,r_name:=paste(r_name,name,sep = '.')][,parent_column:=v2])
      
    }))
    cols_ext[,parent_fill:=color]
    cols_ext[,grandparent_column:=grandparentname]
    
    return(cols_ext)
  }))
  cols=rbind(cols[expanded==0],cols_extended)
  


  return(cols)
}

cols<-fread('analyses_summary/columns_metadata.tsv')[(keep==1)]
colsmtd<-fread('analyses_summary/columns_metametadata.tsv')
colorsmtd<-fread('confidence_colors.tsv')


cols<-PrepColsMtd(cols,colsmtd,res_adxubf)
cols

#get the main sheet
wb<-CreateExcelFormat(res_adxubf,columns_mtd =cols,colors = colorsmtd)


#One sheet per broad context Creation #####
#split per context keeping central information 

res_adxubf_list<-split(res_adxubf,by = 'context_broad2')

for(cont in colsmtd[wildcard=='context_broad2']$r_name){
  message(cont)
  res_adxc<-res_adx[context_broad2==cont]
  
  res_adxcub<-WideTable(res_adxc,group.by = 'qtl_type')
  
    
  #FILTER: keep only variants with GWAS PIP/VCP > 0.1, for maximum of 5
  # if non of the variant has GWAS PIP/VCP > 0.1 we just show top one based on GWAS PIP/VCP and on xQTLPIP/VCP
  
  res_adxcub[,top_variants:=((max_variant_inclusion_probability>=0.1)&(max_variant_inclusion_probability_rank<=5|cV2F_rank<=5))|max_variant_inclusion_probability_rank==1|variant_rank_xqtl==1|cV2F_rank==1,by='ADlocus']
  res_adxcubf<-res_adxcub[(top_variants)]
 
  #add supplemental cols
  res_adxcubf[,gwas_significance:=ifelse(min_pval<5e-8,'genome wide',
                                        ifelse(min_pval<1e-6,
                                               'suggestive',
                                               'ns'))]
  res_adxcubf[,mlog10pval:=-log10(min_pval)]
  
  #variant inclusion top confidence level
  res_adxcubf[,top_confidence:=str_extract(xQTL_effects,'C[0-9]')]
  res_adxcubf$top_confidence|>unique()

  wb<-CreateExcelFormat(res_adxcubf,columns_mtd =cols,colors = colorsmtd,
                    wb = wb,sheet_name = cont)
  
  
}


saveWorkbook(wb, 'Unified_197AD_loci_xQTL_summary_V2.xlsx', overwrite = TRUE)


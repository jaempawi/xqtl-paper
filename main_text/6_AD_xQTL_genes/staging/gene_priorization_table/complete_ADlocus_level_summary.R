 
#Gene locus summary
#need 5 Metadata:
#- metadata_analysis.csv containing all exported tables paths
#- contexts_metadata.csv containing contexts (datasets) label mapping between the different analysis/exported tables
#- columns_metadata.tsv containing all columns to keep for the excel sheets, with associated metadata (column width, coloring..)
#- excel_metadata.tsv containing some meta information for the excel sheet construction
#- pattern_coloring.tsv containing all pattern to color in the excel sheet
#1 metadata is optional
#- long_table_columns_selection.csv to generate a long table with selected columns from the table res_allanalysis_ADloci_overlap.csv.gz generated in step III 
#(in this table each row is a variant-ADlocus-Method-context-gene_name information, and so facilitate querying informations ) 


#install.packages(c('openxlsx'))

source('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/gene_prio_utils.R')

out<-'main_text/6_AD_xQTL_genes/staging/gene_priorization_table/'
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

mtd<-fread('metadata_analysis.csv',header = T)
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

contexts<-fread('contexts_metadata.csv')


res_ts<-fread(file.path('/data',mtd[Method=='trans_finemapping']$Path))

table(res_ts$resource)
cat(sort(unique(res_ts$resource)),sep = '\n')
res_ts<-merge(res_ts,contexts[,.(resource=context_trans,context)],all.x = TRUE)
res_ts[is.na(context)]#ok

table(res_ts$context)
#create locus_context_id for cs95>cs70>cs50
res_ts[,hit:=cs_coverage_0.5!=0|cs_coverage_0.7!=0|cs_coverage_0.95!=0,by=.(event_ID)]

res_tsf<-res_ts[(hit)]
res_tsf[,
        credibleset:=ifelse(cs_coverage_0.95>0,paste0('cs95_',cs_coverage_0.95),
                            ifelse(cs_coverage_0.7>0,paste0('cs70_',cs_coverage_0.7),
                                   paste0('cs50_',cs_coverage_0.5)))]
res_tsf[,region:=gene_ID]

res_tsf[,locuscontext_id:=paste(event_ID,region,credibleset,sep='_')]
res_tsf[,gene_ID:=str_extract(event_ID,'ENSG[0-9]+')]
res_tsf[is.na(gene_ID),.(resource,variant_ID,event_ID,region)]#for gpqtl need mapping
gpmap<-fread('/data/gpQTL/gp_coordinates_clean.txt')
gpmap[str_detect(ensembl_id,';'),ensembl_id:=str_extract(ensembl_id,'ENSG[0-9]+')]
res_tsf[is.na(gene_ID),gp_ID:=str_extract(event_ID,'gp_[0-9]+')]

res_tsf[is.na(gene_ID),gene_ID:=gpmap[gp_ID,on='ID']$ensembl_id]
res_tsf[!is.na(gp_ID)]
res_tsf<-res_tsf[,-'gp_ID']
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



#the sn_sQTL ####

ressqtl<-fread('../../../data/sn-sQTL/major_cell_FDR_0_05.txt.gz',col.names = c('splice_site', 'variant_ID', 'pval', 'beta', 'se', 'FDR', 'if_significant', 'spliced_genes', 'gene_type', 'splice_type', 'cell_type'))
ressqtl
table(ressqtl$splice_type)
table(ressqtl$cell_type)

ressqtl[,gene_name:=strsplit(spliced_genes,',')[[1]][1],by=.(spliced_genes)]
unique(ressqtl$gene_name)
length(unique(ressqtl[is.na(gene_name)]$splice_site))
length(unique(ressqtl$splice_site)) #6895/33512 without gene identified for the splice site
#rm them
ressqtlf<-ressqtl[!is.na(gene_name)]

#get gene id
trans<-fread('/data/resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list')
ressqtlf<-merge(ressqtlf,unique(trans[,.(gene_ID=gene_id,gene_name)]),all.x = TRUE)
ressqtlf[is.na(gene_ID)]$gene_name|>unique() #121/ 
ressqtlf$gene_name|>unique()|>length() #121/ 7238 genes without matching gene name on our reference

ressqtlf<-ressqtlf[!is.na(gene_ID)]

#add harmonized context name
contexts<-fread('contexts_metadata.csv')

ressqtlf<-merge(ressqtlf,unique(contexts[,.(cell_type=context_snsQTL ,context,context_hi)]),by='cell_type')

ressqtlf<-ressqtlf[,-'cell_type']

fwrite(ressqtlf,fp(out,'res_sn_sqtl_summ.csv.gz'))
ressqtlf<-fread(fp(out,'res_sn_sqtl_summ.csv.gz'))

mtd[Method=='LR'&`Data Type`=='sQTL',summary_file:=fp(out,'res_sn_sqtl_summ.csv.gz')]

fwrite(mtd,'metadata_analysis.csv')



# 3)the colocs##### 
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
table(res_cad2$event_ID)
setdiff(unique(res_cad2$event_ID),contexts$context)
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
unique(res_cfs$ge)

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
#first run gwas_cs_extension.R
if(!file.exists(fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))){
  #export table: RUN export_mr.R
  stop('need run gwas_cs_extension.R first to continue')
}

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


#    - MR####
if(!any(mtd$Method=='MR')){
  #export table: RUN export_mr.R
  stop('need run export_mr.R first to continue')
}

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

#cTWAS####
#$cs column not left blank (having values such as L1 L2 L3... 
#although the $pip can be small) and also check $pip column values being greater than 0.75
contexts<-fread(fp(out,'contexts_metadata.csv'))
res_ctwas<-fread(file.path('/data',mtd[Method=='cTWAS']$Path))

res_ctwas[cs!=''&susie_pip>0.75]$molecular_id|>unique() #76
res_ctwas[cs!=''&susie_pip>0.75]$context|>unique()

setdiff(res_ctwas$context,contexts$context)#OK
setnames(res_ctwas,'molecular_id','gene_ID')

res_ctwasf<-res_ctwas[cs!=''&susie_pip>0.75]
res_ctwasf<-res_ctwasf[,-c('group','type')]
res_ctwasf[,cTWAS_signif:=TRUE]
res_ctwasf[(cTWAS_signif)]

fwrite(res_ctwasf,fp(out,'res_AD_cTWAS_pip075.csv.gz'))
res_ctwasf<-fread(fp(out,'res_AD_cTWAS_pip075.csv.gz'))

mtd[Method%in%c('cTWAS'),summary_file:=fp(out,'res_AD_cTWAS_pip075.csv.gz')]


#Metadata of harmonized analysis check and saving####
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

#ok except TWAS/MR for variant_ID and locuscontext_id but normal
mtd[,variant_level_method:=!Method%in%c('MR','twas','cTWAS')]
fwrite(mtd,'metadata_analysis.csv')
mtd<-fread('metadata_analysis.csv')
mtd[summary_file=='']

#II) get unified 'AD associated' loci ####
#outputs:ADlocus variant ADlocus_event ADmethod
#here we want to integrate GWAS locus found in i) single gwas finemapping, ii) ADxQTL coloc, iii) ADxAD colocs

res_gwf<-fread(mtd[Method=='AD_GWAS_finemapping']$summary_file[1])
unique(mtd$Modality)
res_cad<-fread(mtd[Modality=='AD_xQTL_colocalization']$summary_file)[str_detect(context,'^AD')]

res_cad2<-fread(mtd[Modality=='AD_meta_colocalization']$summary_file)
res_cad3<-fread(mtd[Method=='Coloc']$summary_file[1])


res_ad<-rbindlist(list(res_gwf[,.(variant_ID,locuscontext_id,AD_method='single_gwas_finemapping',gwas_source=event_ID,credibleset,PIP)],
                       res_cad[,.(variant_ID,locuscontext_id,AD_method='ADxQTL_coloc',gwas_source=gwas_source,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad2[,.(variant_ID,locuscontext_id,AD_method='ADxAD_coloc',gwas_source=event_ID,vcp,cos_npc,min_npc_outcome,coloc_coef=coef)],
                       res_cad3[,.(variant_ID,locuscontext_id,AD_method='ADFsusie_coloc',gwas_source=gwas_source,SNP_PPH4,L_PP.H4.abf,xQTL_L,AD_L)]),fill = T)
res_ad[,chr:=seqid(variant_ID)]


res_ad<-UnifyLoci(res_ad,variant_col = 'variant_ID',locus_col = 'locuscontext_id',
                  group.by = 'chr')
setnames(res_ad,'uni.locus_id','ADlocus')


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
fwrite(res_ad,fp(out,'AD_loci_unified.csv.gz'))
res_ad<-fread(fp(out,'AD_loci_unified.csv.gz'))


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
res_ad[order(variant_ID,-variant_inclusion_probability),gwas_sources_effects:=paste(gwas_source_effect[!is.na(gwas_source_effect)&!duplicated(str_extract(gwas_source_effect,'[BWJK]'))],collapse =  '|'),by='variant_ID']


# remove variants with min_pval>1e-5 or in cs70 or cs50 only
res_adf<-res_ad[min_pval<1e-5]
unique(res_adf$ADlocusID) #203/249
#only cs95 or coloc
res_adf<-res_adf[ADlocusID%in%ADlocusID[susie_coverage=='cs95'|str_detect(GWAS_methods,'coloc')]]
length(unique(res_adf$ADlocusID)) #197/249
unique(res_adf[min_pval<5e-8]$ADlocusID)|>length() #103
unique(res_adf[susie_coverage=='cs95']$ADlocusID)|>length()  #88
res_adf[is.na(is.cs95),is.cs95:=FALSE]


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
res_adfv<-fread(fp(out,'/AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'))
update_summary_ad=TRUE

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
mismatches<-unique(mismatches)
unique(mismatches$variant_ID)#213/5318 variants mismatch
unique(mismatches[,.(variant_id,variant_ID)])

effect_cols<-c('conditional_effect','coef','effect_direction','top_effect','beta')

mismatches[,to_flip:=ref(variant_id)==alt(variant_ID)|alt(variant_id)==ref(variant_ID)]
mismatches[(to_flip)]
table(mismatches[(to_flip)]$context)
table(mismatches[(to_flip)]$Method)
table(mismatches[(to_flip)][,paste(Method,context)])

#Annotate each variant level Methods table  with the AD Loci
mtd[file.exists(summary_file),summary_file_ad:=paste0(tools::file_path_sans_ext(summary_file[1],compression = T),'_overlapADloci.csv.gz'),by='summary_file']

mtd[file.exists(summary_file)&(variant_level_method),{
    message(Method[1])
    
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
          variants_to_flip<-intersect(res$variant_id,mismatches[(to_flip)]$variant_id)
          
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
      
      
      fwrite(resf,summary_file_ad[1])
    }
    
    summary_file_ad[1]
  },by=c('summary_file')]


fwrite(mtd,'metadata_analysis.csv')
mtd<-fread('metadata_analysis.csv')
res_adfv<-fread(fp(out,'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'),tmpdir = '/adpelle1/tmp/')


#get the long AD overlap table
res_adx<-rbindlist(lapply(unique(mtd[file.exists(summary_file_ad)&(variant_level_method)]$summary_file_ad),function(f){
  message(f)
  m=mtd[summary_file_ad==f]$Method[1]
  message(m)
  if(m=='ColocBoost'){
    m=mtd[summary_file_ad==f]$Modality[1]
    summary_file_ad=mtd[Modality==m]$summary_file_ad[1]
    
  }else if(m=='LR'){
    m='sn_sQTL'
    summary_file_ad=mtd[Method=='LR'&`Data Type`=='sQTL']$summary_file_ad[1]
    
  }else {
    summary_file_ad=mtd[Method==m]$summary_file_ad[1]
    
  }
  res<-fread(summary_file_ad,tmpdir = '/adpelle1/tmp/')[,Method:=m]
  message('found associated with ',length(unique(res$ADlocus)),' AD locus')
  return(res)
}),fill = T)
res_adx
res_adx[!str_detect(gene_ID,'^ENSG')]$context|>table() #OK

res_adx[context=='']$Method|>table() #Coloc and Meta but normal
res_adx[context=='']$context|>table() #Coloc and Meta but normal

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

res_adx[,only_by_proxi:=!any(context%in%c('AD_Bellenguez_EADB_2022','AD_Bellenguez_EADI_2022',
                                          'AD_Wightman_ExcludingUKBand23andME_2021',
                                          'AD_Kunkle_Stage1_2019')),by='locus_index']
res_adx[,APOE_region:=chr==19&pos>43905790 &pos<45905791]
unique(res_adx[(APOE_region)]$locus_index)
#trans[gene_name=="APOE"][,.(start-1e6,end+1e6)]

# #add gwaszscore for gwas method
# res_adx<-merge(res_adx,unique(res_ad[,.(variant_ID,context=gwas_source,gwas_zscore)]),all.x = T,by=c('variant_ID','context'))
res_adx[(is.in.ad.locus)]|>nrow()

#bind to TWAS/MR
#for twas, get ADlocus-gene connection to annotate 
summary_file=mtd[Method=='twas']$summary_file[1]
summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file,compression = T),'_overlapADloci.csv.gz')
update_twasad<-FALSE
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


res_twadf<-fread(summary_file_ad,tmpdir = '/adpelle1/tmp/')
res_twadf[(MR_signif)]
res_twadf[is.na(MR_signif),MR_signif:=FALSE]

#bind it
res_adx<-rbind(res_adx,res_twadf[,Method:='TWAS/MR'],fill=T)

#add cTWAS
summary_file=mtd[Method=='cTWAS']$summary_file[1]
summary_file_ad<-paste0(tools::file_path_sans_ext(summary_file,compression = T),'_overlapADloci.csv.gz')
update_ctwasad<-FALSE
if(!file.exists(summary_file_ad)|update_ctwasad){
  res_adgenes<-unique(res_adx[,.(ADlocus,gene_ID,ADlocusID,locus_index)])
  
  
  res_ctw<-fread(summary_file)
  res_ctwad<-merge(res_ctw[,-c('ADlocus','ADlocusID','locus_index')],res_adgenes,by='gene_ID',all.x = T,allow.cartesian = T)
  res_ctwadf<-res_ctwad[!is.na(ADlocus)]
  unique(res_ctwadf$gene_ID)
  fwrite(res_ctwadf,summary_file_ad)
  mtd[Method%in%c('cTWAS'),summary_file_ad:=..summary_file_ad]
  fwrite(mtd,'metadata_analysis.csv')
}

# #lacking AD genes are of interest ?
# trans<-fread('/data/resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list')
# trans[setdiff(res_ctwad$gene_ID,res_ctwadf$gene_ID),on='gene_id']$gene_name
# # [1] "EPDR1"      "MAP3K1"     "GTPBP1"     "RBX1"       "ZFYVE21"    "RIPK2"      "PLA2G12A"   "HS3ST3B1"   "UBA2"      
# # [10] "COL5A1"     "RPS15A"     "PLXNC1"     "ABCA8"      "ZCCHC24"    "SMG8"       "TRANK1"     "PIK3CD"     "SLC16A11"  
# # [19] "PACS2"      "GPR141"     "XPNPEP3"    "IGHG2"      "AP000295.1" "IFNAR2"     "EEF1G"   
# msig<-fread('/adpelle1/xqtl-paper/resources/all_CPandGOs_gene_and_genesets.csv.gz')
# resor<-OR3(trans[unique(res_ctwad$gene_ID),on='gene_id']$gene_name,split(msig$gene,msig$pathway),background =unique(msig$gene) )
# resor #yes participate to e.g. GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY (PIK3CD,MAP3K1,IGHG2); GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS (RBX1,RIPK2)
# 
# #does the CS of these genes overlap with our AD loci?
# res_ctwadnf<-res_ctwad[!gene_ID%in%res_ctwadf$gene_ID]
# res_ctwasal<-fread('/data/analysis_result/ctwas/archive/export/summary/ctwas_single_fmprs_merged.tsv.gz')
# res_ctwasalf<-merge(res_ctwasal,unique(res_ctwadnf[,.(region_id,cs)]))
# res_ctwasalf[type=='SNP',variant_ID:=paste0('chr',id)]
# unique(res_ctwasalf,by=c('region_id','cs'))
# res_ctwasalfadloci<-merge(res_ctwasalf,res_adfv,by='variant_ID')
# res_ctwadnf[cs=='L1'&region_id=='17_57489969_60570445']
# trans['ENSG00000167447',on='gene_id']#SMG8
# res_adx[locus_index==144&cV2F_rank==1][,.(Method,gene_ID,context)]
# 


#bind it
res_ctwadf<-fread(summary_file_ad,tmpdir = '/adpelle1/tmp/')

res_adx<-rbind(res_adx,res_ctwadf[,Method:='cTWAS'],fill=T)

#populate TWAS info at variant level
res_adx[,TWAS_signif:=any(TWAS_signif),by=.(ADlocus,gene_ID,context,gwas_source)]
res_adx[,MR_signif:=any(MR_signif),by=.(ADlocus,gene_ID,context,gwas_source)]

res_adx[,cTWAS_signif:=any(Method=='cTWAS'),by=.(ADlocus,gene_ID,context)]


table(unique(res_adx,by=c('locuscontext_id','ADlocus','gene_ID','context'))$Method)
# AD_GWAS_finemapping     AD_meta_colocalization     AD_xQTL_colocalization           APOE interaction 
# 341                        249                       2458                         39 
# Coloc         fSuSiE_finemapping           msex interaction  multi_context_finemapping 
# 93                        137                          2                        747 
# multi_gene_finemapping                         QR single_context_finemapping          trans_finemapping 
# 527                        848                       2333                          2 
# TWAS/MR 
# 661 
#add lacking contexts
res_adx[is.na(context)]$locuscontext_id


#split contexts FOR sQTL: if 'UP' put as u-sQTL (Unproductive sQTL), if 'PR' as p-sQTL (productive)
res_adx[str_detect(event_ID,':UP:')&!str_detect(context,'u_sQTL'),context:=str_replace(context,'sQTL','u_sQTL')]
res_adx[str_detect(event_ID,':PR:')&!str_detect(context,'p_sQTL'),context:=str_replace(context,'sQTL','p_sQTL')]
res_adx[str_detect(event_ID,':UP:')]$context


#add the context broad and short
contexts<-fread('contexts_metadata.csv')
setdiff(res_adx$context,contexts$context) 
setdiff(contexts$context,res_adx$context)

res_adx<-merge(res_adx[,-c('context_broad','context_short')],
               unique(contexts[,.(context,context_broad,context_short)]),all.x = T,by='context')
table(res_adx[,.(Method,context)])

table(res_adx[,.(Method,context_short)])


#add GENE info tss /tes
update_geneinfo=TRUE
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
  trans<-fread(fp(out,'genes_infos.csv.gz'),tmpdir = '/adpelle1/tmp/')
}

res_adx<-merge(res_adx[,-c('tss','tes','gene_name')],trans[,-c('start','end','#chr')],all.x = T,by='gene_ID')

res_adx[,chr:=seqid(variant_ID[1],only_num = T),by='variant_ID']
res_adx[,pos:=pos(variant_ID[1]),by='variant_ID']

res_adx[,effect_allele:=alt(variant_ID[1]),by='variant_ID']
res_adx[,distance_from_tss:=pos-tss]
res_adx[,distance_from_tes:=pos-tes]
res_adx[str_detect(gene_ID,';')]$gene_ID|>unique()
res_adx[str_detect(gene_ID,';')]$Method|>table()
res_adx[str_detect(gene_ID,'ENSG')&is.na(gene_name)]
#associate to a gene each epiQTL: populate gene name for epiQTL if geneQTL overlap the locus
epi_methods<-c('fSuSiE_finemapping','Coloc')
geneqtl_methods<-c('single_context_finemapping',
                'multi_context_finemapping','AD_xQTL_colocalization','TWAS/MR')

res_adx<-rbind(res_adx[!Method%in%epi_methods|is.na(Method)],
               res_adx[Method%in%epi_methods]|>split(by=c('locus_index','context'))|>lapply(function(epi){
                 epivars<-epi$variant_ID
                 genes_overlapping_epi<-res_adx[Method%in%geneqtl_methods][variant_ID%in%epivars]$gene_name|>unique()
                 
                 epi_genes<-rbindlist(lapply(genes_overlapping_epi,function(g)epi[,gene_name:=g]))
                 return(epi_genes)
               })|>rbindlist())
res_adx[,gene_ID:=gene_ID[!is.na(gene_ID)][1],by='gene_name']

res_adx[locus_index==1]
table(res_adx$context_broad)
table(res_adx$context_short)


#add broader context and qtl type
res_adx[,context_broad2:=ifelse(context_broad%in%c('bulk_monocyte_eQTL','bulk_macrophage_eQTL','bulk_microglia_eQTL'),'Immune (bulk)',
                                ifelse(str_detect(context_broad,'bulk_brain'),'Brain',
                                       str_remove(context_broad,'_eQTL|_snATAC|_brain|_sQTL')))]
table(res_adx$context_broad2)

res_adx[,qtl_type:=str_extract(context_broad,'(u_|p_)?[a-z]+QTL')]

res_adx[str_detect(context_broad,'snATAC'),qtl_type:='caQTL']
table(res_adx$qtl_type)

#Fix some numerical value
res_adx[twas_z==Inf,twas_z:=max(abs(res_adx[!is.infinite(twas_z)&!is.na(twas_z)][['twas_z']]),na.rm = TRUE)]
res_adx[twas_z==-Inf,twas_z:=-max(abs(res_adx[!is.infinite(twas_z)&!is.na(twas_z)][['twas_z']]),na.rm = TRUE)]


#Summarize table creating column at variant gene -main context level
#here we create the confidence score per variant-gene-context short
unique(res_adx[,.(context,context_short)])

res_adx<-SummarizeTable(res_adx,group.by = 'context_short')
colnames(res_adx)

unique(res_adx[,.(context,context_short)])

fwrite(res_adxa,fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))
table(res_adx$Method)



#OPTIONAL
#keep only columns of interest
# #filter cols: keep only AD locus, gene, method, context, variant_id and gwas infos,  and key method specific infos
colnames(res_adx)|>cat(sep='\n')

if(file.exists('long_table_columns_selection.csv')){

  cols<-fread('long_table_columns_selection.csv')$column_name
  setdiff(cols,colnames(res_adx))
  setdiff(colnames(res_adx),cols)
  
  colsf<-intersect(cols,colnames(res_adx))
  res_adxf<-res_adx[,..colsf]
  res_adxf<-res_adxf[order(locus_index)]
  res_adxf
  
  fwrite(res_adxf,fp(out,'res_allanalysis_ADloci_overlap_selectedcols.csv.gz'))

}





#IV) WIDE TABLE CREATION  ####
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))

res_adxub<-WideTable(res_adx,split.by=c('context_broad2','qtl_type'))

fwrite(res_adxub,fp(out,'res_AD_variants_xQTL.csv.gz'))
res_adxub<-fread(fp(out,'res_AD_variants_xQTL.csv.gz'))

#FILTER: keep only variants with GWAS PIP/VCP > 0.1, for maximum of 5
# if non of the variant has GWAS PIP/VCP > 0.1 we just show top one based on GWAS PIP/VCP and on xQTLPIP/VCP
res_adxub[,top_variants:=((max_variant_inclusion_probability>=0.1)&(max_variant_inclusion_probability_rank<=5|cV2F_rank<=5))|max_variant_inclusion_probability_rank==1|variant_rank_xqtl==1|cV2F_rank==1|(rank(pval)<=1&!is.na(pval)),by='ADlocus']
res_adxubf<-res_adxub[(top_variants)][!is.na(locus_index)][variant_ID!='']
nrow(res_adxubf)#2985
unique(res_adxubf$locus_index)

#some stats
unique(res_adxubf$gene_name)|>length()#1128
unique(res_adxubf[min_pval<5e-8]$gene_name)#602
res_adxubf$top_confidence|>table()
# C1  C2  C3  C4  C5  C6 
# 47  23 197 527 607  11 
#with Trans:
# C1   C2   C3   C4   C5   C6 
# 47   23  197 1366 1443   11 

unique(res_adxubf[order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()

unique(res_adxubf[gwas_significance=='ns'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()


unique(res_adxubf[gwas_significance!='ns'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()
unique(res_adxubf[gwas_significance=='genome wide'][order(locus_index,top_confidence)],
       by='locus_index')$top_confidence|>table()

res_loc<-unique(res_adxubf[!is.na(locus_index)][order(locus_index,top_confidence)],
       by='locus_index')
res_loc[,gwas_significance:=ifelse(gwas_significance=='ns','p<1e-5',ifelse(gwas_significance=='suggestive','p<1e-6','p<5e-8'))]
ggplot(res_loc)+geom_bar(aes(x=top_confidence,fill=gwas_significance))+theme_bw()

res_loc[gwas_significance=='p<5e-8']|>nrow()

#Main Excel Sheet creation ####
#get the columns metadata ready

cols<-fread(fp(out,'columns_metadata.tsv'))[(keep==1)]
colsmtd<-fread(fp(out,'excel_metadata.tsv'))
colorsmtd<-fread(fp(out,'pattern_coloring.tsv'))

cols<-PrepColsMtd(cols,colsmtd,res_adxubf)
unique(cols[,.(parent_column,grandparent_column)])|>tail(100)
unique(cols[,.(parent_column,grandparent_column)])|>tail(100)

setdiff(colnames(res_adxubf),cols$r_name)
#get the main sheet
wb<-CreateExcelFormat(res_adxubf,columns_mtd =cols,
                      wb = wb,colors = colorsmtd)


#One sheet per broad context Creation #####
#split per context keeping central information 

res_adxubf_list<-split(res_adxubf,by = 'context_broad2')

for(cont in colsmtd[wildcard=='context_broad2']$r_name){
  message(cont)
  res_adxc<-res_adx[context_broad2==cont]
  
  
  res_adxc<-SummarizeTable(res_adxc,group.by = 'qtl_type')
  
  res_adxcub<-WideTable(res_adxc)
  
    
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


saveWorkbook(wb, '~/projects-tcwlab/xqtl-resources/data/genes/unified_AD_loci_xQTL_summary.xlsx', overwrite = TRUE)




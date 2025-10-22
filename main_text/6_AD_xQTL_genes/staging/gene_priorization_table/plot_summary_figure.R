out<-'main_text/6_AD_xQTL_genes/staging/gene_priorization_table/'
source('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/gene_prio_utils.R')

#install.packages('ggtext')
 library(ggtext)

#update table figure
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))

res_adx<-res_adx[Method!='trans_finemapping']
#why certain dont have locus id?
res_adx[is.na(locus_index)]$Method|>table()
res_adx$Method|>table()

#add the confidence level

res_adx[,top_confidence:=str_extract(xQTL_effects,'CL[0-9]')]


#xQTL evidence summary dot plot at locus level
#for each QTL evidence type put in 3 cat : absent, present, high, 
#for fp: high if cs95, coloc : high if npc > 0.95, TWAS: MR
table(res_adx$Method)
table(res_adx$context_short)
table(res_adx$context)

#group contexts better: 


res_adx[,locus_gene:=paste(locus_index,ifelse(gene_name=='','?',gene_name),sep='_')]


#annot
res_adx[,genomewide_sig_gene:=any(min_pval<5e-8,na.rm = T),by='gene_name']

res_adx[,n_gwas_gene:=strsplit(gwas_sources,'\\|')|>unlist()|>unique()|>length(),by=.(gene_name)]
res_adx[,n_gwas_locus:=strsplit(gwas_sources,'\\|')|>unlist()|>unique()|>length(),by=.(locus_index)]


#get at locus level and cis information only
res_adxloc<-unique(res_adx[Method!='trans_finemapping'][!context_short%in%c('B','J','K','W')][order(locus_gene,confidence_lvl,cV2F_rank)],by=c('locus_gene','context'))

#top gene per locus
res_adxloc[gene_name!='',n.study.gene.locus:=length(unique(context)),by=.(gene_name,locus_index)]

res_adxloc[gene_name!='',top_gene:=gene_name==gene_name[order(confidence_lvl,-n.study.gene.locus,cV2F_rank)][1],by=.(locus_index)]
res_adxloc[gene_name!='',gene_name[order(confidence_lvl,-n.study.gene.locus,cV2F_rank)][1],by=.(locus_index)]



res_adxlocge<-res_adxloc[context_short!=''][order(locus_index,-cV2F_rank)][(top_gene)]
unique(res_adxlocge$locus_index)|>length() #186
unique(res_adxlocge$gene_name)|>length()#167

#
# group some contexts

res_adxlocge[,context_group:=ifelse(context_short%in%c('bMono eQTL','bMac eQTL','bMic eQTL'),'Immune eQTL',
                                     ifelse(context_short%in%c( 'bulk p-sQTL','bulk u-sQTL','bulk a-sQTL'),'bulk sQTL',
                                            ifelse(context_short%in%c('bulk pQTL','bulk gpQTL'),'bulk (g)pQTL',
                                                   ifelse(context_short%in%c('bulk haQTL','bulk mQTL'),'bulk epiQTL',as.character(context_short)))))]


unique(res_adxlocge$context_group)

unique(res_adxlocge[order(confidence_lvl)],by='gene_name')$confidence_lvl|>table()


# summarize using confidence score  and nstudy
res_adxlocge[,n_study_group:=length(unique(context)),by=.(gene_name,context_group)]
res_adxlocge[,confidence_lvl_group:=sort(confidence_lvl)[1],by=.(gene_name,context_group)]
table(res_adxlocge$confidence_lvl_group)
res_adxlocge[,n_loci_gene:=unique(locus_index)|>length(),by=.(gene_name)]
res_adxlocge[,chr:=seqid(variant_ID[!is.na(variant_ID)][1],only_num = T),by=.(locus_index)]
res_adxlocge[is.na(chr),chr:=`#chr`]
table(res_adxlocge$confidence_lvl_group)


res_adxlocge[,n_study_group.locus:=length(unique(context)),by=.(gene_name,locus_index,context_group)]
res_adxlocge[,confidence_lvl_group.locus:=sort(confidence_lvl)[1],by=.(gene_name,locus_index,context_group)]

res_adxlocge[is.na(APOE_region),APOE_region:=FALSE]

res_adxlocge[,APOE_region.gene:=any(APOE_region),by='gene_name']


fwrite(res_adxlocge,fp(out,'res_summ_top_gene_by_locus.csv.gz'))

res_adxlocge<-fread(fp(out,'res_summ_top_gene_by_locus.csv.gz'))

#genome wide signifcant htis####

contexts_order<-c('Exc eQTL','Exc sQTL',
                  'Inh eQTL','Inh sQTL',
                  'Oli eQTL','Oli sQTL',
                  'OPC eQTL','OPC sQTL',
                  'Ast eQTL','Ast caQTL','Ast sQTL',
                  'Mic eQTL','Mic caQTL',
                  'Immune eQTL',
                  'bulk eQTL',
                  'bulk sQTL',
                  'bulk (g)pQTL',
                  'bulk epiQTL')
setdiff(res_adxlocge$context_group,contexts_order)
res_adxlocge[,context_group:=factor(context_group,levels = contexts_order)]


res_adxlocgef<-res_adxlocge[(genomewide_sig_gene)]
#subset gene with at least CL5
res_adxlocge_cont_top<-unique(res_adxlocgef,by=c('gene_name','context_group'))[gene_name%in%gene_name[confidence_lvl_group%in%c('CL1','CL2','CL3','CL4','CL5')]]

res_adxlocge_cont_top[,gene_name:=factor(gene_name,levels =unique(gene_name[order(chr,pos)]))]

conf_colors<-fread(fp(out,'pattern_coloring.tsv'))[pattern%in%paste0('CL',1:6)]
p<-ggplot(res_adxlocge_cont_top)+
  geom_point(aes(y=gene_name,x=context_group,
                 size=n_study_group,
                 col=confidence_lvl_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+scale_size(range = c(1.5,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = conf_colors$fill_color)
p


#sep causal vs correlated and remove all C6
#C6 is dropped, C1/C2/C3 are called causal, and C4/C5 are called “correlative” or “correlated”
res_adxlocge_cont_topf<-res_adxlocge_cont_top[confidence_lvl_group!='CL6']

#add if frpm, APOE region, or from AD by proxy only
res_adxlocge_cont_topf[,gene_name_2:=as.character(gene_name)]

res_adxlocge_cont_topf[,gene_name_2:=ifelse(any(APOE_region),
                                             paste0("<i>",gene_name_2[1],'</i>'),
                                                    gene_name_2),
                        by='gene_name']

res_adxlocge_cont_topf[,gene_name_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",gene_name_2[1],'</b>'),
                                                    gene_name_2[1]),
                        by='gene_name']

res_adxlocge_cont_topf$gene_name_2|>unique()


res_adxlocge_cont_topf[,gene_name_3:=paste0(gene_name_2,' (',n_gwas_gene,')'),
                        by='gene_name']

res_adxlocge_cont_topf$gene_name_3|>unique()

#add number of loci
res_adxlocge_cont_topf[,gene_name_4:=paste0(gene_name_2,ifelse(n_loci_gene>1,'*',''),' (',n_gwas_gene,')'),
                        by='gene_name']

res_adxlocge_cont_topf$gene_name_4|>unique()

res_adxlocge_cont_topf[,gene_name_4:=factor(gene_name_4,levels = unique(gene_name_4[order(gene_name)]))]

res_adxlocge_cont_topf[is.na(context_group)][,.(context)]|>unique()

#sep in 3categories
res_adxlocge_cont_topf[,confidence_cat_group:=ifelse(confidence_lvl_group%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group=='CL3','Putative causal (CL3)','Associated'))]
res_adxlocge_cont_topf[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated'))]

p<-ggplot(res_adxlocge_cont_topf[!(APOE_region.gene)])+
  geom_point(aes(y=gene_name_4,x=context_group,
                 size=n_study_group,
                 col=confidence_cat_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+
  scale_size(range = c(1.25,5.5),breaks = c(2,7,12))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+  theme(
    axis.text.y = element_markdown(),axis.title.y =element_blank() 
  )
p
ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories.pdf'),height = 9,width = 7)
ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories.png'),height = 9,width = 7)



#version separating each locus####

res_adxlocge_cont_top<-unique(res_adxlocgef,by=c('gene_name','context_group','locus_index'))[gene_name%in%gene_name[confidence_lvl_group.locus%in%c('CL1','CL2','CL3','CL4','CL5')]]

#res_adxlocge_cont_top[,locus_gene:=factor(locus_gene,levels =unique(locus_gene[order(chr,pos)]))]

#sep causal vs correlated and remove all C6
#C6 is dropped
res_adxlocge_cont_topf<-res_adxlocge_cont_top[confidence_lvl_group.locus!='CL6']
#sep in 3categories
res_adxlocge_cont_topf[,confidence_cat_group:=ifelse(confidence_lvl_group.locus%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group.locus=='CL3','Putative causal (CL3)','Associated'))]
res_adxlocge_cont_topf[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated'))]

res_adxlocge_cont_topf[,locus_gene_2:=paste(gene_name,match(locus_gene,
                                                            unique(locus_gene[order(chr,tss,min_pval)])),sep='_'),
                       by=.(gene_name)]


res_adxlocge_cont_topf[,locus_gene_2:=paste0("<i>",locus_gene_2[1],'</i>'),
                       by='locus_gene']

#add number of gwas in this locus and if from genome wide sig loc
res_adxlocge_cont_topf[,locus_gene_2:=paste0(locus_gene_2,' (',n_gwas_locus,ifelse(!genomewide_sig_loc,'*',''),')'),
                       by='locus_gene']


#add if from  AD by proxy only
res_adxlocge_cont_topf[,locus_gene_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",locus_gene_2[1],'</b>'),
                                             locus_gene_2[1]),
                       by='locus_gene']


#order per gene tss and locus position
res_adxlocge_cont_topf[,locus_gene_3:=factor(locus_gene_2,levels = unique(locus_gene_2[order(-chr,-tss,-as.numeric(str_extract(locus_gene_2,'_[0-9+]')|>str_remove('_')))]))]




p<-ggplot(res_adxlocge_cont_topf[!(APOE_region.gene)])+
  geom_point(aes(y=locus_gene_3,x=context_group,
                 size=n_study_group.locus,
                 col=confidence_cat_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+
  scale_size(range = c(1.5,5.5),breaks = c(2,7,12))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+  theme(
    axis.text.y = element_markdown(),axis.title.y =element_blank() 
  )+
  scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
p



ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories.pdf'),height = 12,width = 7)

ggsave(fp(out,'ADloci_xQTL_summary_genome_wide_signif_top_gene_per_locus_3categories.png'),height = 12,width = 7)


#suggestive signal####
#CL1 CL2 CL3 CL4 + more than 2 GWAS + at least 2 datasets in a similar context####


res_adxlocgef<-res_adxlocge[!(genomewide_sig_gene)][gene_name%in%gene_name[confidence_lvl_group%in%c('CL1','CL2','CL3','CL4')&n_gwas_locus>2&n_study_group.locus>=2]]

res_adxlocge_cont_top2f<-unique(res_adxlocgef,by=c('gene_name','context_group','locus_index'))


#sep causal vs correlated 
#add if frpm, APOE region, from AD by proxy only, or from more than one locus

res_adxlocge_cont_top2f[,confidence_cat_group:=ifelse(confidence_lvl_group.locus%in%c('CL1',"CL2"),'Putative causal (CL1, CL2)',ifelse(confidence_lvl_group.locus=='CL3','Putative causal (CL3)','Associated'))]
res_adxlocge_cont_top2f[,confidence_cat_group:=factor(confidence_cat_group,levels = c('Putative causal (CL1, CL2)','Putative causal (CL3)','Associated'))]

#add if frpm, APOE region, or from AD by proxy only
res_adxlocge_cont_top2f[order(locus_index),locus_gene_2:=paste(gene_name,match(locus_index,
                                                                              unique(locus_index[order(chr,tss,pos)])),sep='_'),
                       by=.(gene_name)]


res_adxlocge_cont_top2f[,locus_gene_2:=paste0("<i>",locus_gene_2[1],'</i>'),
                       by='locus_gene']

res_adxlocge_cont_top2f[,locus_gene_2:=ifelse(!all(only_by_proxi,na.rm = T),paste0("<b>",locus_gene_2[1],'</b>'),
                                             locus_gene_2[1]),
                       by='locus_gene']


res_adxlocge_cont_top2f[,locus_gene_3:=paste0(locus_gene_2,' (',n_gwas_locus,')'),
                       by='locus_gene']

res_adxlocge_cont_top2f[,locus_gene_3:=factor(locus_gene_3,levels = unique(locus_gene_3[order(-tss,-pos)]))]



p<-ggplot(res_adxlocge_cont_top2f[!(APOE_region.gene)])+
  geom_point(aes(y=locus_gene_3,x=context_group,
                 size=n_study_group.locus,
                 col=confidence_cat_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+
  scale_size(range = c(1.5,5.5),breaks = c(2,7,12))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = c('brown1','deepskyblue4','darkseagreen3'))+  theme(
    axis.text.y = element_markdown(),axis.title.y =element_blank() 
  )+
  scale_y_discrete(labels=function(x) str_remove(x, "[A-Za-z0-9]+_[2-9]|_[0-9]+"))
p

ggsave(fp(out,'ADloci_xQTL_summary_suggestive_loci_top_gene_per_locus_3categories.png'),height = 5,width = 6)


setdiff(res_adxlocge_cont_top2f$gene_name,res_adxlocge_cont_topf$gene_name)
#v "CCNT2"   "ICA1L"   "BLNK"    "PLEKHA1" "RITA1"   "CTSH" 




setwd('/adpelle1/xqtl-paper-final/')
out<-'main_text/6_AD_xQTL_genes/staging/gene_priorization_table/'
source('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/gene_prio_utils.R')

#install.packages('ggtext')
 library(ggtext)

#update table figure
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))

#add the confidence level

res_adx[,top_confidence:=str_extract(xQTL_effects,'C[0-9]')]


#xQTL evidence summary dot plot at locus level
#for each QTL evidence type put in 3 cat : absent, present, high, 
#for fp: high if cs95, coloc : high if npc > 0.95, TWAS: MR
table(res_adx$Method)
table(res_adx$context_short)
table(res_adx$context)

#group contexts better: 

res_adx[Method%in%'TWAS/MR',evidence_type:='twas']
res_adx[str_detect(Method,'finemapping')&!Method%in%c('AD_GWAS_finemapping','trans_finemapping'),evidence_type:='finemapping']
res_adx[str_detect(str_to_lower(Method),'coloc'),evidence_type:='coloc']

res_adx[,evidence_level:=NULL]
res_adx[evidence_type=='twas',evidence_level:=ifelse(any(MR_signif),2,1),by=.(ADlocus,gene_ID,context)]
res_adx[evidence_type=='finemapping',evidence_level:=ifelse(any(susie_coverage=='cs95'),2,1),by=.(ADlocus,gene_ID,context)]
res_adx[evidence_type=='coloc',evidence_level:=ifelse(any(cos_npc>0.95)|any(L_PP.H4.abf>0.95),2,1),by=.(ADlocus,gene_ID,context)]

res_adx[,locus_gene:=paste(locus_index,ifelse(gene_name=='','?',gene_name),sep='_')]
res_adx[order(locus_gene,cV2F_rank),locus_gene_variant:=paste(locus_gene,variant_ID[1],sep='_'),by='locus_gene']



#Genome wide signif hits
res_adx[,genomewide_sig_loc:=any(min_pval<5e-8),by='locus_gene']


res_adx[context_short=='Ast caQTL']
#get at locus level
res_adxloc<-unique(res_adx[!context_short%in%c('B','J','K','W')][!is.na(evidence_type)][order(locus_gene,confidence_lvl,cV2F_rank)],by=c('context','locus_gene','evidence_type'))

#count numnber of dataset per broad context
res_adxloc[,n.study:=length(unique(context[evidence_level>0])),
        by=c('context_short','locus_gene','evidence_type')]


#Genome wide signif hits

res_adxlocf<-res_adxloc[(genomewide_sig_loc)][evidence_level>0]
unique(res_adxlocf$locus_gene_variant)#241
unique(res_adxlocf$locus_index)|>length() #79
unique(res_adxlocf$Method)

#take top gene per locus
#based on the confidence level

#top gene: based on confidence level and then number of studies/dataset
res_adxlocf[gene_name!='',n.study.gene:=length(unique(context[evidence_level>0])),by=.(gene_name,locus_index)]

res_adxlocf[gene_name!='',top_gene:=gene_name==gene_name[order(top_confidence,-n.study.gene)][1],by=.(locus_index)]



res_adxlocfge<-res_adxlocf[context_short!=''][order(locus_index,-cV2F_rank)][(top_gene)]
unique(res_adxlocfge$locus_index)|>length() #94
unique(res_adxlocfge$gene_name)|>length()#82

#plot the total confidence score per method
res_adxlocfge[gene_name!='',tot.evidence_level:=sum(evidence_level,na.rm = TRUE),by=.(locus_index,gene_name,evidence_type,context_short)]


res_adxlocfge[,chr:=seqid(variant_ID[!is.na(variant_ID)][1],only_num = T),by=.(locus_index)]
res_adxlocfge[is.na(chr),chr:=`#chr`]
res_adxlocfge[is.na(chr)]|>nrow()

res_adxlocfge[,locus_gene:=factor(locus_gene,levels = unique(locus_gene[order(locus_index,confidence_lvl,-n.study.gene)]))]

unique(res_adxlocfge$context_short)
contexts_order<-c('Exc eQTL','Inh eQTL','Oli eQTL','OPC eQTL','Ast eQTL','Ast caQTL','Mic eQTL','bMic eQTL','bMac eQTL','bMono eQTL','bulk eQTL',
                  'bulk p-sQTL','bulk u-sQTL','bulk a-sQTL',
                  'bulk pQTL','bulk gpQTL','bulk mQTL','bulk haQTL')

res_adxlocfge[context_short=='bulk sQTL',context_short:='bulk a-sQTL']


res_adxlocfge[,context_short:=factor(context_short,levels = c(contexts_order))]

p<-ggplot(res_adxlocfge)+geom_point(aes(y=locus_gene,x=evidence_type,
                                        size=tot.evidence_level,
                                        col=as.factor(evidence_level),
                                        shape=evidence_level==2&evidence_type=='twas'))+
  facet_grid(chr~context_short,scales = 'free',space = 'free')+scale_size(range = c(1,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  scale_color_manual(values = c('black','red3'))+
  labs(shape='MR significant')
p
ggsave(fp(out,'genome_wide_signif_ADloci_xQTL_summary.pdf'),height = 14,width = 12)

#with colors by context

p<-ggplot(res_adxlocfge[evidence_level!=0][n.study.gene>0])+geom_point(aes(y=locus_gene,x=evidence_type,
                                                                          size=tot.evidence_level,
                                                                          col=context_short,
                                                                          shape=evidence_level==2&evidence_type=='twas'))+
  facet_grid(chr~context_short,scales = 'free',space = 'free')+scale_size(range = c(1,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(shape='MR significant')
p


ggsave(fp(out,'genome_wide_signif_ADloci_xQTL_summary_contextcolors.pdf'),height = 14,width = 12)


#version by gene

res_adxlocfge[,tot.evidence_level.gene:=sum(evidence_level),by=.(gene_name,evidence_type,context_short)]

res_adxlocfge[,gene_name:=factor(gene_name,levels =unique(gene_name[order(locus_index,confidence_lvl,-n.study.gene)]))]

p<-ggplot(unique(res_adxlocfge[evidence_level!=0][n.study.gene>0],by=c('gene_name','locus_index','evidence_type','context_short')))+
  geom_point(aes(y=gene_name,x=evidence_type,
                 size=tot.evidence_level.gene,
                 col=context_short,
                 shape=evidence_level==2&evidence_type=='twas'))+
  facet_grid(chr~context_short,scales = 'free',space = 'free')+scale_size(range = c(1,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(shape='MR significant')
p

ggsave(fp(out,'genome_wide_signif_ADloci_xQTL_summary_contextcolors_per_gene.pdf'),height = 12,width = 10)


#version grouping some contexts

res_adxlocfge[,context_group:=ifelse(context_short%in%c('bMono eQTL','bMac eQTL','bMic eQTL'),'Immune eQTL',
                                     ifelse(context_short%in%c( 'bulk p-sQTL','bulk u-sQTL','bulk a-sQTL'),'bulk sQTL',
                                     ifelse(context_short%in%c('bulk pQTL','bulk gpQTL'),'bulk (g)pQTL',
                                            ifelse(context_short%in%c('bulk haQTL','bulk mQTL'),'bulk epiQTL',as.character(context_short)))))]
contexts_order<-c('Exc eQTL','Inh eQTL','Oli eQTL','OPC eQTL','Ast eQTL','Ast caQTL','Mic eQTL',
                  'Immune eQTL',
                  'bulk eQTL',
                  'bulk sQTL',
                  'bulk (g)pQTL',
                  'bulk epiQTL')

res_adxlocfge[,context_group:=factor(context_group,levels = contexts_order)]
res_adxlocfge[,tot.evidence_level.gene.contextgroup:=sum(evidence_level),by=.(gene_name,evidence_type,context_group)]

p<-ggplot(unique(res_adxlocfge[evidence_level!=0][n.study.gene>0],by=c('gene_name','locus_index','evidence_type','context_group')))+
  geom_point(aes(y=gene_name,x=evidence_type,
                 size=tot.evidence_level.gene.contextgroup,
                 col=context_group,
                 shape=evidence_level==2&evidence_type=='twas'))+
  facet_grid(chr~context_group,scales = 'free',space = 'free')+scale_size(range = c(1,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(shape='MR significant')
p

ggsave(fp(out,'genome_wide_signif_ADloci_xQTL_summary_contextcolors_per_gene_context_group.pdf'),height = 12,width = 8)


#subset gene to top 2 tiers: with C1/2/3/4 confidence
unique(res_adxlocfge[order(confidence_lvl)],by='gene_name')$confidence_lvl|>table()


res_adxlocfge_cont_top<-unique(res_adxlocfge[evidence_level!=0][n.study.gene>0],by=c('gene_name','locus_index','evidence_type','context_group'))[gene_name%in%gene_name[top_confidence%in%c('C1','C2','C3','C4')]]


p<-ggplot(res_adxlocfge_cont_top)+
  geom_point(aes(y=gene_name,x=evidence_type,
                 size=tot.evidence_level.gene.contextgroup,
                 col=context_group,
                 shape=evidence_level==2&evidence_type=='twas'))+
  facet_grid(chr~context_group,scales = 'free',space = 'free')+scale_size(range = c(0.75,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(shape='MR significant',size='Total evidence level')
p
ggsave(fp(out,'genome_wide_signif_ADloci_xQTL_summary_contextcolors_per_gene_context_group_C1toC4.pdf'),height = 10,width = 8)



# summarize using confidence score instead and nstudy####
res_adxlocfge[,n_study_group:=length(unique(context[evidence_level>0])),by=.(gene_name,context_group)]
res_adxlocfge[,confidence_lvl_group:=sort(confidence_lvl)[1],by=.(gene_name,context_group)]
table(res_adxlocfge$confidence_lvl_group)
res_adxlocfge_cont_top<-unique(res_adxlocfge,by=c('gene_name','context_group'))[gene_name%in%gene_name[top_confidence%in%c('C1','C2','C3','C4')]]


conf_colors<-fread(fp(out,'pattern_coloring.tsv'))[pattern%in%paste0('C',1:6)]
p<-ggplot(res_adxlocfge_cont_top)+
  geom_point(aes(y=gene_name,x=context_group,
                 size=n_study_group,
                 col=confidence_lvl_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+scale_size(range = c(1.5,5))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = conf_colors$fill_color)
p


#sep causal vs correlated
#C6 is dropped, C1/C2/C3 are called causal, and C4/C5 are called “correlative” or “correlated”
library(ggtext)
res_adxlocfge_cont_topf<-res_adxlocfge_cont_top[confidence_lvl_group!='C6']
res_adxlocfge_cont_topf[,confidence_cat_group:=ifelse(confidence_lvl_group%in%c('C1',"C2",'C3'),'putative causal','correlated')]

#add if frpm, APOE region, or from AD by proxy only

res_adxlocfge_cont_topf[is.na(APOE_region),APOE_region:=FALSE]
res_adxlocfge_cont_topf[,gene_name_2:=as.character(gene_name)]

res_adxlocfge_cont_topf[,gene_name_2:=ifelse(any(APOE_region),
                                             paste0("<i>",gene_name_2[1],'</i>'),
                                                    gene_name_2),
                        by='gene_name']

res_adxlocfge_cont_topf[,gene_name_2:=ifelse(!all(only_by_proxi),paste0("<b>",gene_name_2[1],'</b>'),
                                                    gene_name_2),
                        by='gene_name']

res_adxlocfge_cont_topf$gene_name_2|>unique()
res_adxlocfge_cont_topf[is.na(gene_name_2)]
res_adxlocfge_cont_topf[,gene_name_2:=factor(gene_name_2,levels = unique(gene_name_2[order(gene_name)]))]




p<-ggplot(res_adxlocfge_cont_topf)+
  geom_point(aes(y=gene_name_2,x=context_group,
                 size=n_study_group,
                 col=confidence_cat_group))+
  facet_grid(chr~'',scales = 'free',space = 'free')+
  scale_size(range = c(1.25,5.5),breaks = c(2,7,12))+theme_minimal()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  theme(strip.text.x = element_text(angle = 90))+
  labs(size='# datasets',col='Confidence level')+
  scale_color_manual(values = c('cyan3','red'))+  theme(
    axis.text.y = element_markdown()
  )
p






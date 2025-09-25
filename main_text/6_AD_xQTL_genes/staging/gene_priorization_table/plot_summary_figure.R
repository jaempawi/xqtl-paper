setwd('/adpelle1/xqtl-paper-final/')
out<-'main_text/6_AD_xQTL_genes/staging/gene_priorization_table/'
source('../alexandre-utils/r_utils.R')


#update table figure
res_adx<-fread(fp(out,'res_allanalysis_ADloci_overlap.csv.gz'))

#xQTL evidence summary dot plot at locus level
#for each QTL evidence type put in 3 cat : absent, present, high, 
#for fp: high if cs95, coloc : high if npc > 0.95, TWAS: MR
table(res_adx$Method)
table(res_adx$context_short)
#group contexts better: 

res_adx[Method%in%'TWAS/MR',evidence_type:='twas']
res_adx[str_detect(Method,'finemapping')&Method!='AD_GWAS_finemapping',evidence_type:='finemapping']
res_adx[str_detect(str_to_lower(Method),'coloc'),evidence_type:='coloc']

res_adx[,evidence_level:=NULL]
res_adx[evidence_type=='twas',evidence_level:=ifelse(any(MR_signif),2,1),by=.(ADlocus,gene_ID,context_short)]
res_adx[evidence_type=='finemapping',evidence_level:=ifelse(any(susie_coverage=='cs95'),2,1),by=.(ADlocus,gene_ID,context_short)]
res_adx[evidence_type=='coloc',evidence_level:=ifelse(any(cos_npc>0.95)|any(L_PP.H4.abf>0.95),2,1),by=.(ADlocus,gene_ID,context_short)]

res_adx[,locus_gene:=paste(locus_index,ifelse(gene_name=='','?',gene_name),sep='_')]
res_adx[order(locus_gene,cV2F_rank),locus_gene_variant:=paste(locus_gene,variant_ID[1],sep='_'),by='locus_gene']

res_adx[,n.study:=length(unique(context[evidence_level>0])),
            by=c('context_short','locus_gene','evidence_type')]

#Genome wide signif hits
res_adx[,genomewide_sig_loc:=any(min_pval<5e-8),by='locus_gene']

res_adx[context_short=='Ast caQTL']
#get at locus level
res_adxloc<-unique(res_adx[!context_short%in%c('B','J','K','W')][!is.na(evidence_type)][order(locus_gene,cV2F_rank)],by=c('context','locus_gene','evidence_type'))


#Genome wide signif hits

res_adxlocf<-res_adxloc[(genomewide_sig_loc)]
unique(res_adxlocf$locus_gene_variant)#778
unique(res_adxlocf$locus_index)|>length() #98

res_adxlocf[gene_name!='',tot.evidence_level:=sum(evidence_level),by=.(locus_index,gene_name,evidence_type,context_short)]

#take top gene per locus
#top gene: based on n xQTL evidence, or being the in top1 cV2F score
res_adxlocf[gene_name!='',n.xQTL.gene:=sum(n.study),by=.(gene_name)]

res_adxlocf[gene_name!='',top_gene:=gene_name==gene_name[order(-n.xQTL.gene/cV2F_rank)][1],by=.(locus_index)]

res_adxlocfge<-unique(res_adxlocf[context_short!=''][order(locus_index,-cV2F_rank)][(top_gene)])
unique(res_adxlocfge$locus_index)|>length() #95
unique(res_adxlocfge$gene_name)|>length()#83


res_adxlocfge[,chr:=seqid(variant_ID[!is.na(variant_ID)][1],only_num = T),by=.(locus_index,gene_name)]
res_adxlocfge[,pos:=pos(variant_ID)]

res_adxlocfge[,gene_name:=factor(gene_name,levels = unique(gene_name[order(chr,pos)]))]
res_adxlocfge[,locus_gene:=factor(locus_gene,levels = unique(locus_gene[order(chr,pos)]))]

unique(res_adxlocfge$context_short)
contexts_order<-c('Exc eQTL','Inh eQTL','Oli eQTL','OPC eQTL','Ast eQTL','Ast caQTL','Mic eQTL','bMic eQTL','bMac eQTL','bMono eQTL','bulk eQTL',
                  'bulk p-sQTL','bulk u-sQTL','bulk a-sQTL',
                  'bulk pQTL','bulk gpQTL','bulk mQTL','bulk haQTL')
res_adxlocfge[context_short=='bulk sQTL',context_short:='bulk a-sQTL']

setdiff(res_adxlocfge$context_short,contexts_order)
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



#APOE CS
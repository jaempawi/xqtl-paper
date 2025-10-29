#merge_APOE_cs
setwd('/adpelle1/xqtl-paper-final/')
out<-'main_text/6_AD_xQTL_genes/staging/APOElocus'
dir.create(out)
source('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/gene_prio_utils.R')

trans<-fread('/data/resource/references/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.region_list')

resf<-fread('main_text/6_AD_xQTL_genes/staging/gene_priorization_table/res_all_single_context_finemapping_cs50orgreater.csv.gz',tmpdir = '/adpelle1/tmp')

#the 2 merged block###
#chr19:42346101-46842901
resblocks<-resf[chr=='19'&pos>=42346101&pos<=46842901]

resblocks<-merge(resblocks[,gene_id:=gene_ID],trans[,.(gene_id,gene_name)],all.x = T)

#overlap between locus
resblocks[,overlapping_locus.var:=paste(unique(locuscontext_id),collapse ='|'),by='variant_id']

resblocks[,overlapping_locus:=paste(setdiff(unlist(strsplit(overlapping_locus.var,'\\|')),locuscontext_id),collapse = '|'),by='locuscontext_id']
resblocks[,n.overlapping_locus:=length(setdiff(unlist(strsplit(overlapping_locus.var,'\\|')),locuscontext_id)),by='locuscontext_id']

#unify loci
resblocks<-UnifyLoci(resblocks) 

fwrite(resblocks,fp(out,'xqtl_in_APOE_2blocks_region.csv.gz'))

resblocks<-fread(fp(out,'xqtl_in_APOE_2blocks_region.csv.gz'))

cs3<-resblocks[variant_ID=='chr19:44945208:T:G']$uni.locus_id[1]

cs4<-resblocks[variant_ID=='chr19:44946027:T:G']$uni.locus_id[1]

cs3cs4vars<-resblocks[uni.locus_id%in%c(cs3,cs4)]$variant_ID|>unique() #123

#compared to previous version, what is lacking?
resapoeold<-fread('/adpelle1/xqtl-paper/APOElocus/xqtl_in_APOE50kbregion.csv.gz')

setdiff(resapoeold[cs_union_name%in%c('cs3','cs4')]$variant_id,cs3cs4vars)#all are here

# focus with +/-50kb of APOE
apoereg<-get_gene_info(gene_name = 'APOE')
start.reg=apoereg$gene_info$gene_start-50000
end.reg=apoereg$gene_info$gene_end+50000
reg=paste('chr19',start.reg,end.reg,sep='-')


resapoe<-resblocks[chr=='19'&pos>start.reg&pos<end.reg]


fwrite(resapoe,fp(out,'xqtl_in_APOE50kbregion.csv.gz'))

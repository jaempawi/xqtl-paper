
#GWAS extended CS creation
setwd('/adpelle1/xqtl-paper-final/main_text/6_AD_xQTL_genes/staging/gene_priorization_table/')
out<-'.'
source('gene_prio_utils.R')

library(pecotmr)

#if r > 0.8, AND the union of the two CS have min abs(r) > 0.5
#first get long table variant1 variant2 r for abs(r)>0.5

res_gwf <- fread(fp(out,'res_all_single_gwas_finemapping_cs50orgreater.csv.gz'))
res_gwf[,source:='AD_GWAS_finemapping']
#bind with ADxQTL, ADAD and fsusie colocs
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))
res_cad<-fread(fp(out,'res_coloc_meta_AD.csv.gz'))
res_cad[,chr:=seqid(variant_ID)]
res_cad[,pos:=pos(variant_ID)]

res_cfsf<-fread(fp(out,'res_coloc_AD_epiQTL.csv.gz'))
res_cfsf[,.(variant_ID,cos_ID,context)]
res_cfsf[,chr:=seqid(variant_ID)]
res_cfsf[,pos:=pos(variant_ID)]

res_cosad<-rbindlist(list(res_c[str_detect(context,'^AD')][,source:='AD_xQTL_colocalization'],
                 res_cad[,source:='AD_meta_colocalization'],
                 res_cfsf[,source:='AD_epiQTL_colocalization']),fill=TRUE)
all_variants<-Reduce(union,list(res_gwf$variant_ID,res_c$variant_ID,res_cad$variant_ID,res_cfsf$variant_ID))
sum(duplicated(seqid(all_variants))&duplicated(pos(all_variants)))#ok

#for 1 block
# b1<-res_gwf$region[1]
# ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
#                         region = data.frame(chr=seqid(b1),
#                                             start=start(b1),
#                                             end=end(b1)))$combined_LD_matrix
# variants<-intersect(rownames(ldblock),res_gwf$variant_id)
# variants_dt<-data.table(variant_id=variants)
# variants_dt
# variants_cors<-variants_dt[,{
#   vars_cor<-rownames(ldblock)[abs(ldblock[,variant_id])>0.5]
#   
#   list(variant2=vars_cor,
#        r=ldblock[vars_cor,variant_id])
# },by=.(variant_id)]
# variants_cors

#for all block
ldmeta<-fread('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv')

regions=bed_inter(data.table(chr=seqid(all_variants),
                             start=pos(all_variants)-1,
                             end=pos(all_variants),
                             variant_ID=all_variants),
                  ldmeta[,.(`#chrom`,start,end,block=paste(`#chrom`,start,end,sep='_'))])[[8]]|>unique()


variants_cors<-rbindlist(lapply(1:length(regions),function(i){
  message(i,'/',length(regions))
  b=regions[i]
  ldblock<-load_LD_matrix('/data/resource/ADSP_R4_EUR_LD/ld_meta_file.tsv',
                          region = data.frame(chr=seqid(b),
                                              start=start(b),
                                              end=end(b)))$combined_LD_matrix
  variants<-intersect(rownames(ldblock),str_remove(all_variants,'chr'))
  variants_dt<-data.table(variant_id=variants)
  variants_cors<-variants_dt[,{
    vars_cor<-rownames(ldblock)[abs(ldblock[,variant_id])>0.5]
    
    list(variant2=vars_cor,
         r=ldblock[vars_cor,variant_id])
  },by=.(variant_id)]
  variants_cors[,region:=b][]
}))
variants_cors[,variant_ID:=paste0('chr',variant_id)]
variants_cors[,variant_ID2:=paste0('chr',variant2)]

fwrite(variants_cors,fp(out,'gwas_variants_cor0.5.csv.gz'))
variants_cors<-fread(fp(out,'gwas_variants_cor0.5.csv.gz'))

#if r > 0.8, AND the union of the two CS have min abs(r) > 0.5
#per region, iterate in all CS of a GWAS finding if correlation with other GWAS CS meet this condition, if it meet,
#extend the CS 
variants_corsf<-variants_cors[variant_ID%in%all_variants&variant_ID2%in%all_variants]
#iteration 1
i<-1
gwas_new<-rbindlist(lapply(regions,function(r){
  message(r)
  res_gwff<-res_gwf[region==r]
  #get long table of the CS
  res_gwff_cs<-melt(res_gwff,measure.vars=c('cs_coverage_0.95',
                                                'cs_coverage_0.7',
                                                'cs_coverage_0.5'),
                         variable.name='coverage',value.name = 'cs_num')
  res_gwff_cs<-res_gwff_cs[cs_num!=0]
  
  #add the coloc sets
  res_cf<-res_cosad[chr==seqid(r)&pos>=start(r)&pos<=end(r)]

  res_gwff_cs<-rbind(res_gwff_cs,res_cf,fill=TRUE)

  #iterate per CS

  res_gwff_cs_ext<-rbindlist(lapply(unique(res_gwff_cs$locuscontext_id),function(l){
    message(l)
    res_gwff_cs_study<-res_gwff_cs[locuscontext_id==l]
    res_gwff_cs_others<-res_gwff_cs[locuscontext_id!=l]
    res_gwff_cs_study_ext<-res_gwff_cs_study[,{
      vars_cs<-variant_ID #the CS
      #check if this CS have any cor >0.8 with CS from others GWAS and mincorr 0.8
      #if yes, extract these variants
      to_add<-res_gwff_cs_others[,to_merge:={
        other_cs<-variant_ID
        mincorr50_vars<-variants_corsf[variant_ID%in%vars_cs&variant_ID2%in%other_cs]
        if(nrow(mincorr50_vars)>0){
          #any 0.8 between those CS
          have_highLD<-mincorr50_vars[,any(abs(r)>0.8)]
          #all variants pairs have r>0.5
          have_all50<-nrow(mincorr50_vars)==length(other_cs)*length(vars_cs)
          rep(have_highLD&have_all50,.N) #rep is to return TRUE or FALSE for each Variants of thise CS
        }else{
          rep(FALSE,.N)
        }
       
      },
      by=.(locuscontext_id)][(to_merge)]$variant_ID
      
      new<-setdiff(to_add,vars_cs)
      if(length(new)>0){
        message(length(new),' variants added')

      }
      
      merged_cs<-union(vars_cs,to_add)
      
      data.table(variant_ID=merged_cs)[,original_cs:=variant_ID%in%vars_cs]
      
    },by=.(locuscontext_id,source)]
    
      return(res_gwff_cs_study_ext)
    }),fill = TRUE)
    return(res_gwff_cs_ext[,region:=r])
  }),fill = TRUE)


#iteration 2 and +
change<-TRUE
while(change){
  i<-i+1
  message('iteration ', i)
  change=FALSE
  
  gwas_toploci_old<-copy(gwas_new)
  gwas_new<-rbindlist(lapply(regions,function(r){
    message(r)
    res_gwff_cs<-gwas_new[region==r]
   
    #iterate by study
    res_gwff_cs_ext<-rbindlist(lapply(unique(res_gwff_cs$locuscontext_id),function(l){
      message(l)
      res_gwff_cs_study<-res_gwff_cs[locuscontext_id==l]
      res_gwff_cs_others<-res_gwff_cs[locuscontext_id!=l]
      #iterate per CS
      res_gwff_cs_study_ext<-res_gwff_cs_study[,{
        vars_cs<-variant_ID #the CS
        vars_original<-variant_ID[(original_cs)]
        #check if this CS have any cor >0.8 with CS from others GWAS and mincorr 0.8
        #if yes, extract these variants
        to_add<-res_gwff_cs_others[,to_merge:={
          other_cs<-variant_ID
          mincorr50_vars<-variants_corsf[variant_ID%in%vars_cs&variant_ID2%in%other_cs]
          if(nrow(mincorr50_vars)>0){
            #any 0.8 between those CS
            have_highLD<-mincorr50_vars[,any(abs(r)>0.8)]
            #all variants pairs have r>0.5
            have_all50<-nrow(mincorr50_vars)==length(other_cs)*length(vars_cs)
            rep(have_highLD&have_all50,.N) #rep is to return TRUE or FALSE for each Variants of thise CS
          }else{
            rep(FALSE,.N)
          }
          
        },
        by=.(locuscontext_id)][(to_merge)]$variant_ID
        new<-setdiff(to_add,vars_cs)
        if(length(new)>0){
          message(length(new),' variants added')
          change<<-TRUE
          
        }
        merged_cs<-union(vars_cs,to_add)
        
        data.table(variant_ID=merged_cs)[,original_cs:=variant_ID%in%vars_original]
        
      },by=.(locuscontext_id,source)]
      return(res_gwff_cs_study_ext)
    }),fill = TRUE)
    return(res_gwff_cs_ext[,region:=r])
  }),fill = TRUE)

  #gwas_new[,cs_num_ext:=ifelse(gwas_toploci_old,cs_num,-cs_num)]
  ntotadded=nrow(gwas_new[!(original_cs)])
  message(ntotadded,' total variants added')

}


message('converged after ',i,' iteration')
gwas_new
#recreate top_loci table
fwrite(gwas_new,fp(out,'all_adlocis_extended_any0.8ANDmin0.5union_withfsusie.csv.gz'))
gwas_new<-fread(fp(out,'all_adlocis_extended_any0.8ANDmin0.5union_withfsusie.csv.gz'))

#saved back the different methods
#FINEMAP
table(gwas_new$source)
gwas_newf<-gwas_new[source=='AD_GWAS_finemapping']
gwas_newf[,study:=strsplit(locuscontext_id,'_chr')[[1]][1],by='locuscontext_id']
gwas_newf[,coverage:=ifelse(str_detect(locuscontext_id,'cs95'),'cs_coverage_0.95',
                            ifelse(str_detect(locuscontext_id,'cs70'),'cs_coverage_0.7',
                                   'cs_coverage_0.5'))]
gwas_newf[,cs_num:=as.numeric(str_extract(locuscontext_id,'[0-9]+$'))]

resfp<-dcast(gwas_newf,
                study+region+variant_ID~coverage,
                value.var = 'cs_num')
#add pip, z
mtd<-fread('metadata_analysis.csv',header = T,select = 1:6)

res_gw<-rbindlist(lapply(file.path('/data',mtd[Method=='AD_GWAS_finemapping']$Path),function(f)fread(f)),fill = T)
res_gw[,study:=event_ID]
resfp<-merge(resfp,res_gw[,.(variant_ID,PIP,z,study,cs_coverage_0.95_original=cs_coverage_0.95,
                                   cs_coverage_0.7_original=cs_coverage_0.7, cs_coverage_0.5_original=cs_coverage_0.5)],
                by=c('study','variant_ID'),all.x = T)
resfp[,event_ID:=study]
resfp[is.na(cs_coverage_0.95),cs_coverage_0.95:=0]
resfp[is.na(cs_coverage_0.7),cs_coverage_0.7:=0]
resfp[is.na(cs_coverage_0.5),cs_coverage_0.5:=0]
resfp[is.na(cs_coverage_0.95_original),cs_coverage_0.95_original:=0]
resfp[is.na(cs_coverage_0.7_original),cs_coverage_0.7_original:=0]
resfp[is.na(cs_coverage_0.5_original),cs_coverage_0.5_original:=0]

#distrib n.extension 
resfp[,n.extension:=sum(cs_coverage_0.95!=cs_coverage_0.95_original),by=.(study,region,cs_coverage_0.95)]
summary(unique(resfp[cs_coverage_0.95>0],by=c('study','region','cs_coverage_0.95'))$n.extension)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.000   0.000   2.000   7.559   6.000 144.000 

resfp[,n.extension:=sum(cs_coverage_0.7!=cs_coverage_0.7_original),by=.(study,region,cs_coverage_0.7)]
summary(unique(resfp[cs_coverage_0.7>0],by=c('study','region','cs_coverage_0.95'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.00    0.00    5.00   18.99   23.00  144.00 

resfp[!is.na(PIP)]

#export
#fwrite(resfp,'/adpelle1/export/AD_GWAS_finemapping_109_blocks_top_loci_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz')

fwrite(resfp,fp(out,'res_all_single_gwas_finemapping_cs50orgreater_unified_withAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#COLOCS
table(gwas_new$source)

#ADmeta coloc
res_cad<-fread(fp(out,'res_coloc_meta_AD.csv.gz'))
res_cad[,chr:=seqid(variant_ID)]
res_cad[,pos:=pos(variant_ID)]


res_cadnew<-gwas_new[source=='AD_meta_colocalization']

res_cadnew<-merge(res_cadnew,res_cad,all.x = T,
                  by=c('locuscontext_id','variant_ID'))
res_cadnew[,chr:=seqid(variant_ID)]
res_cadnew[,pos:=pos(variant_ID)]

res_cadnew<-rbindlist(lapply(unique(res_cad$event_ID),function(g){
  loci<-unique(res_cosad[event_ID==g]$locuscontext_id)
  
  rbind(res_cadnew[locuscontext_id%in%loci&!is.na(event_ID)],res_cadnew[locuscontext_id%in%loci&is.na(event_ID),event_ID:=g][])
  
}))
res_cadnew[,context:=event_ID]
res_cadnew[,gwas_source:=event_ID]
res_cadnew<-unique(res_cadnew)
res_cadnew[is.na(locuscontext_id)]
table(res_cadnew$locuscontext_id)
unique(res_cadnew$locuscontext_id)

#extension
res_cadnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
summary(unique(res_cadnew,by=c('locuscontext_id'))$n.extension)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   0.0     0.0     2.0    21.4    11.0  1036.0 


fwrite(res_cadnew,fp(out,'res_coloc_meta_AD_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#ADxQTL coloc
res_c<-fread(fp(out,'res_coloc_AD_xQTL_unified.csv.gz'))

table(gwas_new$source)
table(gwas_new$source)

res_cnew<-gwas_new[source=='AD_xQTL_colocalization']
unique(res_cnew$locuscontext_id)
res_cnew[,context:=paste0('AD_',strsplit(locuscontext_id,'AD_|_ENSG')[[1]][2]),by='locuscontext_id']
unique(res_cnew$context)

res_cnew<-merge(res_cnew,res_c,all = T,
                  by=c('locuscontext_id','context','variant_ID'))

# res_cnew[,event_ID:=event_ID[!is.na(event_ID)][1],by=.(locuscontext_id,context)]
res_cnew[,context_coloc:=context_coloc[!is.na(context_coloc)][1],by=.(locuscontext_id,context)]
#res_cnew[,context:=context[!is.na(context)][1],by=.(locuscontext_id,context)]

res_cnew[,gwas_source:=gwas_source[!is.na(gwas_source)][1],by=.(locuscontext_id)]


#extension
res_cnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
summary(unique(res_cnew,by=c('locuscontext_id'))$n.extension)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #  0.0     3.0    10.0   135.2    49.0  2639.0 

res_cnew[is.na(locuscontext_id)]
fwrite(res_cnew,fp(out,'res_coloc_AD_xQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))

#fsusie Coloc
res_cfsf<-fread(fp(out,'res_coloc_AD_epiQTL.csv.gz'))
table(res_cfsf$context)

table(gwas_new$source)

res_cfnew<-gwas_new[source=='AD_epiQTL_colocalization']
res_cfnew<-merge(res_cfnew,unique(res_cfsf[,.(locuscontext_id,gwas_source)]),by='locuscontext_id')

res_cfnew<-merge(res_cfnew,res_cfsf,all = T,
                by=c('locuscontext_id','gwas_source','variant_ID'))

# res_cfnew[is.na(context)]
# res_cfnew[,context:=context[!is.na(context)][1],by=.(locuscontext_id)]

#extension
res_cfnew[,n.extension:=sum(is.na(cos_ID)),by='locuscontext_id']
summary(unique(res_cfnew,by=c('locuscontext_id'))$n.extension)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 0.00    0.00    3.00   28.31   36.00  284.00 

res_cfnew[is.na(locuscontext_id)]
fwrite(res_cfnew,fp(out,'res_coloc_AD_epiQTL_unified_withFP_andAllCoS_any0.8ANDmin0.5_converged.csv.gz'))



#QC####
# #how many GWAS overlap compared to before?
# gwas_new[,overlapping_locus.var:=paste(unique(locuscontext_id),collapse ='|'),
#          by=c('variant_ID')]
# 
# 
# gwas_new[,overlapping_locus:=paste(sort(unique(unlist(strsplit(overlapping_locus.var,'\\|')))),collapse = '|'),
#          by=c('locuscontext_id')]
# gwas_new[,n.overlapping_locus:=str_count(overlapping_locus,'AD'),
#          by=c('locuscontext_id')]
# gwas_new[n.overlapping_locus>0]$locuscontext_id|>unique()|>length()#606/859
# unique(gwas_new$locuscontext_id)|>length()
# 
# res_gwf_cs<-melt(res_gwf,measure.vars=c('cs_coverage_0.95_min_corr',
#                                       'cs_coverage_0.7_min_corr',
#                                       'cs_coverage_0.5_min_corr'),
#                 variable.name='coverage',value.name = 'cs_num')
# res_gwf_cs<-unique(res_gwf_cs[cs_num!=0])
# res_gwf_cs[,locuscontext_id:=paste(region,study,coverage,cs_num,sep='_')]
# res_gwf_cs[,overlapping_locus.var:=paste(unique(locuscontext_id),collapse ='|'),
#          by=c('variant_ID')]
# 
# res_gwf_cs[,overlapping_locus:=paste(sort(unique(unlist(strsplit(overlapping_locus.var,'\\|')))),collapse = '|'),
#          by=c('locuscontext_id')]
# res_gwf_cs[,n.overlapping_locus:=str_count(overlapping_locus,'AD'),
#          by=c('locuscontext_id')]
# res_gwf_cs[n.overlapping_locus>0]$locuscontext_id|>unique()|>length()#575/859
# unique(res_gwf_cs$locuscontext_id)|>length()#86
# 
# sum(unique(res_gwf_cs,by='locuscontext_id')$n.overlapping_locus)#2719
# 
# sum(unique(gwas_new,by='locuscontext_id')$n.overlapping_locus)#3045
# 
# #cs95 only
# res_gwf_loc<-unique(res_gwf_cs[coverage=='cs_coverage_0.95_min_corr'],by='locuscontext_id')
# 
# res_gwf_loc[n.overlapping_locus>1]$locuscontext_id|>unique()|>length()#117
# 
# res_gwf_loc$locuscontext_id|>unique()|>length()#117/170
# 
# 
# gwas_new_loc<-unique(gwas_new[coverage=='cs_coverage_0.95_min_corr'],by='locuscontext_id')
# 
# gwas_new_loc[n.overlapping_locus>1]$locuscontext_id|>unique()|>length()#117
# 
# gwas_new_loc$locuscontext_id|>unique()|>length()#117/170
# 
# #total overlap for cs95
# sum(res_gwf_loc$n.overlapping_locus)-nrow(res_gwf_loc)#373
# sum(gwas_new_loc$n.overlapping_locus)-nrow(gwas_new_loc)#402
# 
# #how many loci compared to before?
# gwas_new_loc[coverage=='cs_coverage_0.95_min_corr']$overlapping_locus|>unique()|>length()#92
# gwas_new_loc
# 
# #before [from]
# res_gwf_loc[coverage=='cs_coverage_0.95_min_corr']$overlapping_locus|>unique()|>length()#95
# res_gwf_cs[,chr:=paste0('chr',seqid(variant_id)),by='locuscontext_id']
# 
# res_gwf_cs[,gwas_region.hit:=paste(chr,min(pos(variant_id)),max(pos(variant_id)),sep='_'),by='locuscontext_id']
# res_gwf_cs[,gwas_region.hit.start:=start(gwas_region.hit)]
# res_gwf_cs[,gwas_region.hit.end:=end(gwas_region.hit)+1]
# 
# fwrite(unique(res_gwf_cs[,.(chr,gwas_region.hit.start,gwas_region.hit.end,locuscontext_id)])[order(chr,gwas_region.hit.start)],'gwas.region.old.bed',sep='\t',col.names=F)
# 
# #intersect regions
# cmds<-list(paste('bedtools merge -i','gwas.region.old.bed',
#                  '>','gwas.region.old.merge.bed'),
#            paste('bedtools intersect -wa -wb -a','gwas.region.old.bed',
#                  '-b', 'gwas.region.old.merge.bed',
#                  '>','gwas.region.old.inter.bed'))
# for(cmd in cmds){
#   system(cmd)
# }
# gwas.inter.bed<-fread('gwas.region.old.inter.bed',select=4:7,col.names = c('locuscontext_id','chr','start','end'))
# gwas.inter.bed[,gwas_region.hit.bedunion:=paste(chr,start,end,sep='-')]
# res_gwf_loc<-unique(res_gwf_cs,by='locuscontext_id')
# 
# 
# res_gwf_loc<-merge(res_gwf_loc[,-'gwas_region.hit.bedunion'],
#                 unique(gwas.inter.bed[,.(locuscontext_id,gwas_region.hit.bedunion)]),by='locuscontext_id')
# #find overlapping AD CS
# res_gwf_loc[,gwas_region.hit.union:={
#   message('======= new locus =======')
#   print(gwas_region.hit.bedunion)
#   #get the different sets list
#   #set number back to the locus
#   #for each gwasoverlap, if intersect, i create the union 
#   overlaps_list<-strsplit(overlapping_locus,'\\|')
#   overlaps_list<-lapply(1:.N, function(i)c(overlaps_list[[i]],locuscontext_id[i]))
#   overlap_set_num<-1:.N
#   already_set<-c()
#   for(i in 1:.N){
#     if(!i%in%already_set){
#       overlap_set_num[i]<-i
#       already_set<-c(already_set,i)
#     }
#     if(.N>1){
#       for(j in (i+1):.N){
#         #if next CSs are intersecting with i, merge with i 
#         if(length(intersect(overlaps_list[[overlap_set_num[i]]],
#                             overlaps_list[[overlap_set_num[j]]]))>0){
#           message('found intersection between ',locuscontext_id[i],' and ',locuscontext_id[j])
#           overlaps_list[[overlap_set_num[i]]]<-union(overlaps_list[[overlap_set_num[i]]],overlaps_list[[overlap_set_num[j]]])
#           
#           overlap_set_num[j]<-overlap_set_num[i]
#           message('new sets :')
#           print(overlap_set_num)
#           already_set<-c(already_set,j)
#           #update the overlap list
#           
#           
#           
#         }
#         
#       }
#     }
#     
#   }
#   paste(gwas_region.hit.bedunion,overlap_set_num,sep='_')
# },by='gwas_region.hit.bedunion']
# 
# res_gwf_loc[str_detect(coverage,'95')]$gwas_region.hit.union|>unique()|>length()
# 
# 
# #n.extension
# ggplot(gwas_new)+geom_bar(aes(x=coverage,fill=original_cs))+theme_bw()
# 
# gwas_new[,n.extension:=.N-sum(original_cs),by='locuscontext_id']
# 
# ggplot(unique(gwas_new,by='locuscontext_id'))+
#   geom_boxplot(aes(x=coverage,y=n.extension))+
#   theme_bw()+theme_bw()



#now add to landscape
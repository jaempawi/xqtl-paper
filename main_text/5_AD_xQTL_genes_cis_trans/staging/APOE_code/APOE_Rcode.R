library(tidyverse)
library(data.table)
library(vroom)
library(vctrs)
library(susieR)
library(matrixStats)
library(colocboost)
library(pecotmr)


APOE_summary <- readRDS("xqtl_only_APOE_all_cohorts_addGWAS.rds")
APOE_summary <- APOE_summary %>% filter(start >= 42346101 & start <= 46842901)

# - summary statistics
# - those are the LD meta path for AD GWAS
studies = c("AD_Bellenguez", "AD_Kunkle_Stage1_2019", "AD_Wightman_Full_2021", "AD_Wightman_Excluding23andMe_2021", 
                    "AD_Wightman_ExcludingUKBand23andME_2021", "AD_Bellenguez_EADB")
sumstat_path_list = c('/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/AD_Bellenguez_2022.sumstats.tsv.gz',
                      '/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/Kunkle_etal_Stage1_results.txt_file_1_hg38.sorted.txt.gz',
                      '/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/PGCALZ2full.hg38.sorted.txt.gz',
                      '/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/PGCALZ2sumstatsExcluding23andMe.hg38.sorted.txt.gz',
                      '/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.hg38.sorted.txt.gz',
                      '/mnt/jast/hpc/gaowang/users/xc2270/colocboost/pipeline/sumstat/EADB_core.tsv.gz')
column_file_path_list = c('/home/xc2270/COLOCBoost/pipeline/Bellenguez_new.yml',
                          '/mnt/vast/hpc/csg/data_public/GWAS_sumstats/GWAS_sumstat_column_mapping/Kunkle_stage_1.yml',
                          '/mnt/vast/hpc/csg/data_public/GWAS_sumstats/GWAS_sumstat_column_mapping/AD_Wightman_Full_2021.yml',
                          '/mnt/vast/hpc/csg/data_public/GWAS_sumstats/GWAS_sumstat_column_mapping/AD_Wightman_Excluding23andMe_2021.yml',
                          '/mnt/vast/hpc/csg/data_public/GWAS_sumstats/GWAS_sumstat_column_mapping/AD_Wightman_ExcludingUKBand23andME_2021.yml',
                          '/mnt/vast/hpc/csg/data_public/GWAS_sumstats/GWAS_sumstat_column_mapping/EADB.yml')
# Replace _regions with _meta_info[1] which is the association window region
n_samples = c('0','0','0','0','0','0') %>% as.numeric
n_cases = c('111326','21982','90338','86531','39918','20301') %>% as.numeric
n_controls = c('677663','41944','1036225','676386','358140','21839') %>% as.numeric

gwas_meta <- data.frame(
    "study_id" = studies,
    "file_path" = sumstat_path_list,
    "column_mapping_file" = column_file_path_list,
    "n_sample" = n_samples,
    "n_case" = n_cases,
    "n_control" = n_controls
)


ldblock<-load_LD_matrix('ld_meta_file_apoe.tsv', region = data.frame(chr='chr19', start=42346101, end=46842901))

impute_opts=list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)
qc_method='rss_qc'
association_window = "chr19:42346101-46842901"
sumstats_region_name_col = NULL
comment_string = NULL
extract_sumstats_region_name = NULL


cos <- APOE_summary$identifier %>% unique
cos %>% length
source("merge_coloc_also_within_loci.R")
threshold <- 0.2
flat_cos <- lapply(1:length(cos), function(i){
    pos <- which(APOE_summary$identifier == cos[i])
    APOE_summary$variant_ID[pos]
})
flat_cos_vcp <- lapply(1:length(cos), function(i){
    pos <- which(APOE_summary$identifier == cos[i])
    APOE_summary$vcp[pos] %>% as.numeric
})
merge_pairwise_idx <- get_merge_pairwise_idx(flat_cos, flat_cos_vcp, threshold = threshold)
flat_coloc_outcome <- lapply(1:length(cos), function(i){
    pos <- which(APOE_summary$identifier == cos[i])
    tmp <- APOE_summary$event_ID[pos] %>% unique
    tmp <- lapply(tmp, function(tt) strsplit(tt, "; ")[[1]] )
    tmp %>% unlist %>% unique
})
    
    
length(merge_pairwise_idx)
    
final_cos <- final_vcp <- list()
final_colocOutcome <- c()
for (ii in 1:length(merge_pairwise_idx)){
    p.merge <- merge_pairwise_idx[[ii]]
    # - coloc outcomes 
    oo <- flat_coloc_outcome[p.merge] %>% unlist
    colocOutcome <- paste0(unique(oo), collapse = "; ")
    # - coloc CoS and pph4
    snps <- unlist(flat_cos[p.merge])
    vcps <- unlist(flat_cos_vcp[p.merge])
    context_df <- data.frame(SNP = snps, vcp = vcps, stringsAsFactors = FALSE)
    unique_snps <- unique(context_df$SNP)
    max_vcp <- sapply(unique_snps, function(snp) {
      max(context_df$vcp[context_df$SNP == snp], na.rm = TRUE)
      # min(context_df$vcp[context_df$SNP == snp], na.rm = TRUE)
    })
    merged_df <- data.frame(SNP = unique_snps, MaxVCP = max_vcp, stringsAsFactors = FALSE)
    cos <- merged_df$SNP
    vcp <- merged_df$MaxVCP
    # - coloc purity
    final_cos <- c(final_cos, list(cos))
    final_vcp <- c(final_vcp, list(vcp))
    final_colocOutcome <- c(final_colocOutcome, colocOutcome)
}
names(final_cos) <- 
    names(final_vcp) <- 
    names(final_colocOutcome) <- 
    paste0("ind_cos_", 1:length(merge_pairwise_idx))
    

gwas_idx <- (_GWAS_)
cos_idx <- (_IDX_)
    
gwas <- gwas_meta$study_id[gwas_idx]
message(gwas)
sumstat_path=gwas_meta %>% filter(study_id==gwas) %>% pull(file_path)
column_file_path=gwas_meta %>% filter(study_id==gwas) %>% pull(column_mapping_file)
n_sample=gwas_meta %>% filter(study_id==gwas) %>% pull(n_sample)
n_case=gwas_meta %>% filter(study_id==gwas) %>% pull(n_case)
n_control=gwas_meta %>% filter(study_id==gwas) %>% pull(n_control)

#get summstats
rss_input <- load_rss_data(
  sumstat_path = sumstat_path, column_file_path = column_file_path,
  n_sample = n_sample, n_case = n_case, n_control = n_control,
  region = association_window
)
message(rss_input$n)
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
qc_results <- summary_stats_qc(preprocess_results$sumstats, list(combined_LD_matrix=preprocess_results$LD_mat), n = rss_input$n,
                             var_y = rss_input$var_y,
                             method = qc_method)


# Perform imputation 
# Leave one CoS out
cos_idx <- final_cos[[cos_idx]]
pp_cos <- which( paste0("chr", preprocess_results$sumstats$variant_id) %in% cos_idx )
    
if (length(pp_cos) != 0){
    
    sumstats_obs <- qc_results$sumstats
    qc_results$sumstats <- qc_results$sumstats[-pp_cos, ]
    impute_results <- raiss(ldblock$ref_panel, qc_results$sumstats, ldblock$combined_LD_matrix,
                  rcond = impute_opts$rcond,
                  R2_threshold = impute_opts$R2_threshold,
                  minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb)
    
    ll <- list(summstats=data.table(impute_results$result_filter),
               sumstats_obs = sumstats_obs,
               outlier_number=qc_results$outlier_number)
    ff <- paste0('results/', gwas, '/rss_qc_imputed_',gwas,'_indcos_', cos_idx, '.rds')

    saveRDS(ll, file.path(ff))

}



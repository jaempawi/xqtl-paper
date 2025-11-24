library(tidyverse)
library(data.table)
library(vroom)
library(vctrs)
library(susieR)
library(matrixStats)
library(pecotmr)
library(colocboost)
merged_info <- readRDS("xqtl_only_APOE_all_cohorts_merged_cos_cs_after_between_purity.rds")
cos <- merged_info$final_set

message("A total of ", cos %>% length, " indepdent sets after checking between purity.")


# - summary statistics
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


ldblock<-load_LD_matrix('ld_meta_file_apoe.tsv',
                         region = data.frame(chr='chr19', start=42346101, end=46842901))

impute_opts=list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)
qc_method='slalom'
association_window = "chr19:42346101-46842901"
sumstats_region_name_col = NULL
comment_string = NULL
extract_sumstats_region_name = NULL


# Leave one CoS out
gwas_idx <- (_GWAS_)
cos_idx <- (_IDX_)


gwas = gwas_meta$study_id[gwas_idx]
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
res_slalom <- slalom(zScore = preprocess_results$sumstats$z,
                     LD_mat = preprocess_results$LD_mat,r2_threshold=0.2)
qc_results<-list()
qc_results$sumstats<-merge(data.table(preprocess_results$sumstats),
                           data.table(res_slalom$data,keep.rownames = 'variant_id')[,.(variant_id,nlog10p_dentist_s)],
                           by='variant_id')
qc_results$LD_mat<-preprocess_results$LD_mat[qc_results$sumstats$variant_id,qc_results$sumstats$variant_id]
qc_results$outlier_number<-nrow(data.table(res_slalom$data,keep.rownames = 'variant_id')[(outliers)])


cs_variant <- cos[[cos_idx]]
message("Current set includs ", length(cs_variant), " variants.")
cs_idx <- gsub("chr", "", cs_variant)
sumstats_obs <- qc_results$sumstats
variants_to_include <- sumstats_obs$variant_id
zscores = sumstats_obs %>% select(chrom, pos, variant_id, A1, A2, z) %>% as.data.frame()
zscores_rm <- zscores %>% filter(!(variant_id %in% cs_idx))
message("Number of variants before removing set is: ", nrow(zscores), ".\n", 
        "Number of variants after removing set is: ", nrow(zscores_rm), ".")

imputed_zscores<-raiss(ref_panel = data.table(ldblock$ref_panel)[variants_to_include,on='variant_id']|>as.data.frame(),
                     known_zscores =zscores_rm,
                     LD_matrix=ldblock$combined_LD_matrix[variants_to_include,variants_to_include],
                     R2_threshold=0.4,minimum_ld = 0.5)

zscores_obs <- zscores %>% filter(variant_id %in% cs_idx)
zscores_rssqc_imputed <- imputed_zscores$result_filter %>% filter(variant_id %in% zscores$variant_id)
summstats_qced <- qc_results$sumstats %>% filter(variant_id %in% zscores$variant_id)
data <- data.frame(
    variant_id = zscores_obs$variant_id,
    z_observed = zscores_obs$z,
    z_imputed = zscores_rssqc_imputed$z[match(zscores_obs$variant_id, zscores_rssqc_imputed$variant_id)],
    z_qced = summstats_qced$z[match(zscores_obs$variant_id, summstats_qced$variant_id)]
)

dir_folder <- paste0("/raise_results/", gwas)
if (!dir.exists(dir_folder)) {
  dir.create(dir_folder, recursive = TRUE)
} 
ff <- paste0(dir_folder, '/rss_qc_imputed_',gwas,'_indcos_', i, '.rds')
saveRDS(data, file.path(ff))
    

get_flat <- function(contexts_info){
    
    flat_info <- list()
    flat_names <- c()
    for (main_name in names(contexts_info)) {
    sublist <- contexts_info[[main_name]]
    if (is.list(sublist)) {
      for (sub_name in 1:length(sublist)) {
        flat_info <- c(flat_info, sublist[sub_name])
        flat_names <- c(flat_names, main_name)
      }
    } else {
      flat_info <- c(flat_info, sublist)
      flat_names <- c(flat_names, main_name)
    }
    }
    names(flat_info) <- flat_names
    return(flat_info)
    
}

get_flat_info <- function(res_contexts){
    
    contexts_snps <- lapply(res_contexts, function(tm){
        snp <- tm$colocalized_variables
        snp_sets <- list()
        for (ii in 1:nrow(tm)){
            ss <- strsplit(snp[ii], "; ")
            snp_sets <- c(snp_sets, ss)
        }
        snp_sets
    })
    contexts_pph4 <- lapply(res_contexts, function(tm){
        pph4 <- tm$colocalized_variables_pph4
        pph4_sets <- list()
        for (ii in 1:nrow(tm)){
            ss <- as.numeric(unlist(strsplit(pph4[ii], "; ")))
            pph4_sets <- c(pph4_sets, list(ss))
        }
        pph4_sets
    })
    contexts_purity <- lapply(res_contexts, function(tm){ sapply(tm$purity, list) })
    contexts_outcomes <- lapply(res_contexts, function(tm){ sapply(tm$colocalized_outcomes, list) })
    
    flat_cos <- get_flat(contexts_snps)
    flat_pph4 <- get_flat(contexts_pph4)
    flat_purity <- get_flat(contexts_purity)
    flat_outcomes <- get_flat(contexts_outcomes)
    
    flat_info <- list("flat_cos" = flat_cos,
                      "flat_pph4" = flat_pph4,
                      "flat_purity" = flat_purity,
                      "flat_outcomes" = flat_outcomes)
    return(flat_info)
    
}


get_merge_pairwise_idx <- function(flat_cos, flat_pph4, threshold = 0.5){
    
    ncsets <- length(flat_cos)
    if (ncsets == 1){ return(lapply(1:ncsets, function(jj) jj)) }
    if_overlap <- matrix(0, nrow =  ncsets, ncol =  ncsets)
    for (i in 1:(ncsets-1)){
        for (j in (i+1):ncsets){
            cset1 <- flat_cos[[i]]
            cset2 <- flat_cos[[j]]
            over_snp <- intersect(cset1, cset2)
            pos1 <- match(over_snp, cset1)
            pos2 <- match(over_snp, cset2)
            pph4_1 <- flat_pph4[[i]][pos1]
            pph4_2 <- flat_pph4[[j]][pos2]
            prop1 <- sum(pph4_1) / sum(flat_pph4[[i]])
            prop2 <- sum(pph4_2) / sum(flat_pph4[[j]])
            # if_overlap[i,j] <- if_overlap[j,i] <- length(intersect(cset1, cset2))!=0
            if_overlap[i,j] <- if_overlap[j,i] <- (prop1 >= threshold) && (prop2 >= threshold)
        }
    }
    
    merge_sets <- function(vec) {
        split_lists <- lapply(vec, function(x) as.numeric(unlist(strsplit(x, ";"))))
        result <- list()
        while (length(split_lists) > 0) {
            current <- split_lists[[1]]
            split_lists <- split_lists[-1]
            repeat {
                overlap_index <- NULL
                for (i in seq_along(split_lists)) {
                    if (length(intersect(current, split_lists[[i]])) > 0) {
                        overlap_index <- i
                        break
                    }
                }
                if (!is.null(overlap_index)) {
                    current <- union(current, split_lists[[overlap_index]])
                    split_lists <- split_lists[-overlap_index]
                } else {
                    break
                }
            }
            result <- c(result, list(paste(sort(current), collapse = ";")))
        }
        return(result)
    }

    if (sum(if_overlap)!=0){
        
        temp <- sapply(1:nrow(if_overlap), function(x){
            tt <- c(x, which(if_overlap[x,] != 0))
            return(paste0(sort(tt), collapse = ";"))
        })
        temp <- merge_sets(temp)
        potential_merged <- lapply(temp, function(x) as.numeric(unlist(strsplit(x, ";"))))
    
    } else {
        potential_merged <- lapply(1:ncsets, function(jj) jj)
    }

  return(potential_merged)
}
                                   
get_union_pairwise <- function(res_contexts, threshold = 0.5) {
    
    flat_info <- get_flat_info(res_contexts)
    flat_cos <- flat_info$flat_cos
    flat_pph4 <- flat_info$flat_pph4
    flat_purity <- flat_info$flat_purity
    flat_outcomes <- flat_info$flat_outcomes
    
    merge_pairwise_idx <- get_merge_pairwise_idx(flat_cos, flat_pph4, threshold = threshold)
    final_cos <- final_pph4 <- list()
    final_colocPhen <- final_purity <- final_colocOutcome <- c()
    for (ii in 1:length(merge_pairwise_idx)){
        p.merge <- merge_pairwise_idx[[ii]]
        # - coloc phenotype
        colocPhen <- paste0(unique(names(flat_cos[p.merge])), collapse = "; ")
        # - coloc outcomes 
        oo <- unlist(lapply(flat_outcomes[p.merge], function(ot) strsplit(ot, "; ")[[1]]))
        colocOutcome <- paste0(unique(oo), collapse = "; ")
        # - coloc CoS and pph4
        snps <- unlist(flat_cos[p.merge])
        pph4s <- unlist(flat_pph4[p.merge])
        context_df <- data.frame(SNP = snps, pph4 = pph4s, stringsAsFactors = FALSE)
        unique_snps <- unique(context_df$SNP)
        max_pph4 <- sapply(unique_snps, function(snp) {
          max(context_df$pph4[context_df$SNP == snp], na.rm = TRUE)
        })
        merged_df <- data.frame(SNP = unique_snps, Maxpph4 = max_pph4, stringsAsFactors = FALSE)
        cos <- merged_df$SNP
        pph4 <- merged_df$Maxpph4
        # - coloc purity
        purity <- min(unlist(flat_purity[p.merge]))
        
        final_cos <- c(final_cos, list(cos))
        final_pph4 <- c(final_pph4, list(pph4))
        final_colocPhen <- c(final_colocPhen, colocPhen)
        final_colocOutcome <- c(final_colocOutcome, colocOutcome)
        final_purity <- c(final_purity, purity)
    }
    names(final_cos) <- final_colocPhen
    names(final_pph4) <- final_colocPhen
    names(final_purity) <- final_colocPhen
    names(final_colocOutcome) <- final_colocPhen
    
    final_info <- list("CoS_union" = final_cos,
                       "CoS_union_outcomes" = final_colocOutcome,
                       "Max_PPH4_union" = final_pph4)
    return(final_info)
    
}
                            
                            
                            

merge_coloc_sets <- function(same_region_path, overlap_pip_sum_threshold = 0.5){
    
    coloc_results <- lapply(same_region_path, readRDS)
    names(coloc_results) <- basename(same_region_path)
    
    ncos <- length(coloc_results)
    coloc_summary <- list()
    for (i in 1:ncos){
        tmp <- coloc_results[[i]][[1]]
        analysis_region <- tmp$analysis_region
        colocalized_outcomes <- paste0(gsub(".rds", "", basename(tmp$analysis_parameter$file1)), "_", tmp$analysis_parameter$context1, 
                                       "; ", 
                                       gsub(".rds", "", basename(tmp$analysis_parameter$file2)), "_", tmp$analysis_parameter$context2)
        purity <- tmp$sets$purity$min.abs.corr
        coloc_idx <- tmp$sets$true_summary %>% rownames %>% as.integer
        colocalized_variables <- sapply(tmp$sets$cs, paste0, collapse = "; ")
        num_variables <- sapply(tmp$sets$cs, length)
        snp_pph4 <- tmp$results %>% select(coloc_idx+1)
        highest_pph4 <- as.numeric(apply(snp_pph4, 2, max))
        colocalized_variables_pph4 <- sapply(1:length(tmp$sets$cs), function(ii){
            snp <- tmp$sets$cs[[ii]]
            pos <- match(snp, tmp$results$snp)
            cos_pph4 <- snp_pph4 %>% select(all_of(ii)) %>% unlist %>% as.numeric
            paste0(cos_pph4[pos], collapse = "; ")
        })
        coloc_summary[[i]] <- data.frame(
            colocalized_outcomes = colocalized_outcomes,
            purity = purity,
            num_variables = num_variables,
            highest_pph4 = highest_pph4,
            colocalized_variables = colocalized_variables,
            colocalized_variables_pph4 = colocalized_variables_pph4
        )
    }
    names(coloc_summary) <- names(coloc_results)
    
    merged_cos = get_union_pairwise(coloc_summary, threshold = overlap_pip_sum_threshold)

    output <- list(
        "coloc_files" = same_region_path,
        "cos_summary" = coloc_summary,
        "merged_cos" = merged_cos
    )
    return(output)
}
#' Gets a list of files containing ASE information for a given tissue
#'
#'
#'
#' @param fullName Tissue name
#' @return dataframe
#' @import tidyverse
#' @export
#'
get_ASE_files_list <- function(fullName) {
  # Read list with samples available
  list <- readr::read_delim("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data/list.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  sp <- function(x) { strsplit(x, '[.]')[[1]][1] }
  list$id <- purrr::map_chr(list$X1, sp)
  colnames(list)[1] <- "file_name"
  # Read data with counts so that we know tissues ids
  path_rda <- "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/data_Rda/"
  pathLoad <- paste(path_rda, fullName, sep="")
  pathLoad <- paste(pathLoad, ".Rda", sep="")
  load(pathLoad)
  ids_tissue <- data.frame(SAMPID=as.character(subset_tissue$SAMPID))
  sp <- function(x) {
    a = unlist(strsplit(x, '[-]')[[1]])
    return(paste(a[1], a[2], sep="-"))
  }
  ids_tissue$id <- purrr::map_chr(as.character(ids_tissue$SAMPID), sp)
  # Merge
  ids_tissue_merged <- merge(ids_tissue, list, by.x="id", by.y="id")
  return(ids_tissue_merged)
}

#' Gets a table with reference and allelic counts for a given sample and given tissue
#'
#'
#'
#' @param sample_id sample ID of interest
#' @param tissue tissue of interest
#' @return dataframe
#' @import tidyverse
#' @export
#
get_counts_for_sample <- function(sample_id, tissue) {
  # Open file with the ASE data for a given sample
  f_name <- paste0("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data/phe000014.v1.GTEx_MidPoint_Imputation_ASE.expression-matrixfmt-ase.c1/", sample_id[1])
  f_name <- paste0(f_name, ".ase_table.tsv.gz")
  df <- readr::read_delim(f_name, "\t", escape_double = FALSE, trim_ws = TRUE)
  # Subset to the tissue of interest
  df <- df %>% filter(TISSUE_ID == tissue) %>% filter(LOW_MAPABILITY ==0 & MAPPING_BIAS_SIM == 0 & GENOTYPE_WARNING ==0)
  # Check if the tissue is present
  if (dim(df)[1] == 0) return(NULL)
  else {
    df <- df[,c(21,1:5,9,10,19,20)]
    df <- df[!is.na(df$GENE_ID),]
    df$sw <- ifelse(df[,9]== "GT;1|0", 1, 0)
    df$REF_COUNT_UPD <- df$REF_COUNT
    df$ALT_COUNT_UPD <- df$ALT_COUNT
    df[df$sw==1,"REF_COUNT_UPD"] <- df[df$sw==1,"ALT_COUNT"]
    df[df$sw==1,"ALT_COUNT_UPD"] <- df[df$sw==1,"REF_COUNT"]
    df <- df[,c(1,12,13)]
    colnames(df)[2] <- paste0("ref_", gsub("-", "\\.", as.character(unlist(sample_id[2]))))
    colnames(df)[3] <- paste0("alt_", gsub("-", "\\.", as.character(unlist(sample_id[2]))))
  }
  return(df)
}

#' Aggregates reads counts from SNPs to transcript and calculates ratios
#'
#' @param df dataframe with count in columns 2 and 3
#' @return dataframe
#
aggregate_counts <- function(df, threshold = 10) {
  df_aggr <- aggregate(df[,c(2,3)], by=list(df$GENE_ID), sum)
  colnames(df_aggr)[1] <- "GENE_ID"
  df_aggr <- df_aggr[df_aggr[,2] + df_aggr[,3] >= threshold,]
  df_aggr$ratio <- df_aggr[,2] / (df_aggr[,2] + df_aggr[,3])
  colnames(df_aggr)[4] <- paste0("ratio", substr(colnames(df_aggr)[2], 4, nchar(colnames(df_aggr)[2])))
  return(df_aggr)
}

#' Reads counts for a given tissue, aggregates them and puts together ratios for all samples
#'
#' @param fullName full tissue name
#' @return dataframe
#
read_counts_and_calc_ratios <- function(fullName) {
  # Get samples ids
  ids_tissue <- get_ASE_files_list(fullName)
  # Get tissue code
  tissue_code <- get_tissue_code(fullName)
  # Read and aggregate the data
  df <- get_counts_for_sample(ids_tissue[1,], tissue_code)
  df_agg <- aggregate_counts(df)
  df_agg <- df_agg[,c(1,4)]
  for (i in 2:dim(ids_tissue)[1]) {
    tmp <- get_counts_for_sample(ids_tissue[i,], tissue_code)
    tmp <- aggregate_counts(tmp)
    if (!is.null(tpm)) df_agg <- merge(df_agg, tmp[,c(1,4)], by.x="GENE_ID", by.y="GENE_ID", all.x = TRUE, all.y = TRUE)
  }
  return(df_agg)
}

#' Reads counts for a given tissue, aggregates them and puts together counts for all samples
#'
#' @param fullName full tissue name
#' @return dataframe
#
read_counts_and_merge_samples <- function(fullName) {
  # Get samples ids
  ids_tissue <- get_ASE_files_list(fullName)
  # Get tissue code
  tissue_code <- get_tissue_code(fullName)
  # Create dataframe for output
  df_agg <- data.frame(GENE_ID = character(), ref_count = numeric(), alt_count = numeric(), sample_id = numeric())
  # Read and aggregate the data
  for (i in 1:dim(ids_tissue)[1]) {
    tmp <- get_counts_for_sample(ids_tissue[i,], tissue_code)
    tmp <- aggregate_counts(tmp)
    tmp <- tmp[,c(1:3)]
    colnames(tmp)[2:3] <- c("ref_count", "alt_count")
    tmp$sample_id <- ids_tissue$SAMPID[i]
    if (!is.null(tpm)) df_agg <- rbind(df_agg, tmp)
  }
  return(df_agg)
}

# To test
#fullName <-  "Kidney - Cortex"
#read_counts_and_calc_ratios(fullName)
#read_counts_and_merge_samples(fullName)

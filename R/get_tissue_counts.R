#' Proccesses GTEx counts and merges counts statistics with classification data
#'
#'
#'
#' @param name Tissue name
#' @param keepCounts Keep columns with counts in the output
#' @param sample_names Default is set to all meaning "all samples"; can be a string of names (c("GTEX-111CU-1326-SM-5NQ8L", "GTEX-111YS-1426-SM-5GID8", "GTEX-1122O-1326-SM-5H11F"))
#' @param mode "norm_counts" / "log" meaning using normalized by library size counts or using log(normalized_counts)
#' @param path_to_class Path to file name_class.txt with MAE/BAE classification call (output by MAGIC or RNA-seq analysis)
#' @param path_rda Path to the folder contatining all .Rda files for tissues
#' @return dataframe with MAE/BAE classification and counts statistics
#' @export
process_to_counts <- function(name,
                              keepCounts = FALSE,
                              sample_names = "all",
                              mode = "norm_counts",
                              path_to_class = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/classification/",
                              path_rda = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/data_Rda/") {
  # Get fullName by name
  fullName <- get_full_name(name)
  if (is.na(fullName)) {
    print("Full name not found, try again")
    return()
  }
  geneSt <- get_classification(name, path_to_class)
  geneSt <- geneSt[!duplicated(geneSt$Ensembl.Gene.ID),]
  # Load counts data
  sv <- load_counts(fullName, path_rda)
  # Subset to samples
  if (sample_names[1]!="all") {
    sv <- sv[, sampleNames(sv) %in% sample_names]
  }
  # Calculate counts statistics
  dataNorm <- get_counts_statistics(sv, keepCounts = keepCounts, mode = mode)
  # Merge classification and statistics
  dataAll <- merge(geneSt, dataNorm, by.x = "Ensembl.Gene.ID", by.y = "geneId")
  return(dataAll)
}

#' Loads counts from Rda
#'
#'
#'
#' @param fullName tissue name
#' @param path_rda Path to the folder contatining all .Rda files for tissues
#' @return ExpressionSet
#' @export
#'
load_counts <- function(fullName,
                        path_rda = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/data_Rda/") {
  tryCatch({
    pathLoad <- paste(path_rda, fullName, sep="")
    pathLoad <- paste(pathLoad, ".Rda", sep="")
    load(pathLoad)
    sv <- subset_tissue
    print(paste0("Number of samples for ", fullName, ":", dim(sv)[2]))
  }, warning = function(war) {
    print(".Rda file not found")
    return()
  })
  return(sv)
}

#' Calculates mean, sd, CV, CV2 for each gene
#'
#'
#'
#' @param sv ExpressionSet for a given tissue
#' @param min_counts_per_gene Minimum number of counts per gene to be considered
#' @param keepCounts Keep columns with counts in the output
#' @param mode "norm_counts" / "log" meaning using normalized by library size counts or using log(normalized_counts)
#' @return dataframe with counts statistics
#' @export
#'
get_counts_statistics <- function(sv, min_counts_per_gene = 10,
                                  keepCounts = FALSE,
                                  mode = "norm_counts") {
  sv_counts <- data.frame(exprs(sv))
  sv_counts <- sv_counts[rowSums(sv_counts)>=min_counts_per_gene, ]
  sv_counts <- sv_counts[rowSums(sv_counts>0)>(dim(sv)[2]/2), ]
  lib.size <- DESeq2::estimateSizeFactorsForMatrix(sv_counts)
  nmat <- as.matrix(t(t(sv_counts)/lib.size))
  if (mode == "log") nmat <- log(nmat)
  dataNorm <- data.frame(nmat)
  dataNorm$rowSd <- matrixStats::rowSds(nmat)
  dataNorm$rowMean <- rowMeans(nmat)
  dataNorm$rowVar <- matrixStats::rowVars(nmat)
  dataNorm$CV <- matrixStats::rowSds(nmat) / rowMeans(nmat)
  dataNorm$CV2 <- dataNorm$CV * dataNorm$CV
  #quant_threshold <- quantile(dataNorm$CV2, 0.99)
  #dataNorm <- dataNorm[dataNorm$CV2<=quant_threshold,]
  dataNorm$gene <- row.names(dataNorm)
  dataNorm$geneId <- matrix(unlist(strsplit(dataNorm$gene,"\\.")), ncol = 2 , byrow = TRUE )[,1]
  if (keepCounts) return(dataNorm)
    else return(data.frame(dataNorm[,c("geneId", "rowMean", "rowSd", "rowVar", "CV", "CV2")]))
}

#' Proccesses GTEx counts and merges counts statistics with classification data
#'
#'
#'
#' @param name Tissue name
#' @param path_to_class Path to file name_class.txt with MAE/BAE classification call (output by MAGIC or RNA-seq analysis)
#' @param path_rda Path to the folder contatining all .Rda files for tissues
#' @return dataframe with MAE/BAE classification and counts statistics
#' @export
process_to_counts <- function(name,
                            path_to_class = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/classification/",
                            path_rda = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/data_Rda/") {
  # Get fullName by name
  fullName <- get_full_name(name)
  if (is.na(fullName)) {
    print("Full name not found, try again")
    return()
  }
  # Read classification, convert to Ensembl ids
  convertionTable <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/convertionTable.txt")
  tryCatch({
    pathToStatus <- paste(path_to_class, name, sep="")
    pathToStatus <- paste(pathToStatus, "_class.txt", sep="")
    genes_status <- read.csv(pathToStatus, sep="\t")
  }, warning = function(war) {
    print("Classification file not found")
    return()
  })
  geneSt <- merge(convertionTable, genes_status, by.x="WikiGene.Name", by.y="Gene")
  # Load counts data
  tryCatch({
    pathLoad <- paste(path_rda, fullName, sep="")
    pathLoad <- paste(pathLoad, ".Rda", sep="")
    load(pathLoad)
    sv <- subset_tissue
    print("Dimentions of the data:")
    print(paste0(dim(sv)[1], " transcripts"))
    print(paste0(dim(sv)[2], " samples"))
  }, warning = function(war) {
    print(".Rda file not found")
    return()
  })
  # Calculate counts statistics
  dataNorm <- get_counts_statistics(sv)
  # Merge classification and statistics
  dataAll <- merge(geneSt, dataNorm, by.x = "Ensembl.Gene.ID", by.y = "geneId")
  return(dataAll)
}

#' Calculates mean, sd, CV, CV2 for each gene
#'
#'
#'
#' @param sv ExpressionSet for a given tissue
#' @param min_counts_per_gene Minimum number of counts per gene to be considered
#' @return dataframe with counts statistics
#' @export
#'
get_counts_statistics <- function(sv, min_counts_per_gene = 10) {
  sv_counts <- data.frame(exprs(sv))
  sv_counts <- sv_counts[rowSums(sv_counts)>=min_counts_per_gene, ]
  sv_counts <- sv_counts[rowSums(sv_counts>0)>(dim(sv)[2]/2), ]
  lib.size <- DESeq2::estimateSizeFactorsForMatrix(sv_counts)
  nmat <- as.matrix(t(t(sv_counts)/lib.size))
  dataNorm <- data.frame(nmat)
  dataNorm$rowSd <- matrixStats::rowSds(nmat)
  dataNorm$rowMean <- rowMeans(nmat)
  dataNorm$rowVar <- matrixStats::rowVars(nmat)
  dataNorm$CV <- matrixStats::rowSds(nmat) / rowMeans(nmat)
  dataNorm$CV2 <- dataNorm$CV * dataNorm$CV
  dataNorm[dataNorm$CV2<=quantile(dataNorm$CV2, 0.99),]
  dataNorm$gene <- row.names(dataNorm)
  dataNorm$geneId <- matrix(unlist(strsplit(dataNorm$gene,"\\.")), ncol = 2 , byrow = TRUE )[,1]
  return (data.frame(dataNorm[,c("geneId", "rowMean", "rowSd", "rowVar", "CV", "CV2")]))
}

#' Gets full tissue name from a short name
#'
#'
#'
#' @param name Short tissue name
#' @return full tissue name
#' @export
#'
get_full_name <- function(name) {
  tissueNames <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/tissueNames.txt", stringsAsFactors = F)
  fullNames <- tissueNames$fullName
  names <- tissueNames$name
  names_list <- setNames(as.list(fullNames), names)
  fullName <- unlist(names_list[name])[1]
  if (is.null(fullName)) return(NA)
  return(fullName)
}

#' Gets tissue code from a short name
#'
#'
#'
#' @param name Short tissue name
#' @return tissue code
#' @export
#'
get_tissue_code <- function(name) {
  tissueNames <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/tissueNames.txt", stringsAsFactors = F)
  fullNames <- tissueNames$fullName
  codes <- tissueNames$code
  codes_list <- setNames(as.list(codes), fullNames)
  code <- unlist(codes_list[name])[1]
  if (is.null(code)) return(NA)
  return(code)
}

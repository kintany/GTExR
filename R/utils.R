#' Gets full tissue name from a short name
#'
#'
#'
#' @param name Short tissue name
#' @return full tissue name
#' @export
#'
get_full_name <- function(name) {
  tissueNames <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/util_files/tissueNames.txt", stringsAsFactors = F)
  fullNames <- tissueNames$fullName
  names <- tissueNames$name
  names_list <- setNames(as.list(fullNames), names)
  fullName <- unlist(names_list[name])[1]
  if (is.null(fullName)) return(NA)
  return(fullName)
}

#' Gets short tissue name from a full name
#'
#'
#'
#' @param fullName Full tissue name
#' @return short tissue name
#' @export
#'
get_short_name <- function(fullName) {
  tissueNames <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/util_files/tissueNames.txt", stringsAsFactors = F)
  fullNames <- tissueNames$fullName
  names <- tissueNames$name
  names_list <- setNames(as.list(names), fullNames)
  name <- unlist(names_list[fullName])[1]
  if (is.null(name)) return(NA)
  return(name)
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
  tissueNames <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/util_files/tissueNames.txt", stringsAsFactors = F)
  fullNames <- tissueNames$fullName
  codes <- tissueNames$code
  codes_list <- setNames(as.list(codes), fullNames)
  code <- unlist(codes_list[name])[1]
  if (is.null(code)) return(NA)
  return(code)
}

#' Returns tissue classification with Ensembl id column
#'
#'
#' @param path_to_class  Path to file name_class.txt with MAE/BAE classification call (output by MAGIC or RNA-seq analysis)
#' @param name Short tissue name
#' @return tissue code
#' @export
#'
get_classification <- function(name, path_to_class = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/classification/") {
  # Read classification, convert to Ensembl ids
  convertionTable <- read.delim("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/convertionTable.txt")
  tryCatch({
    pathToStatus <- paste0(path_to_class, name, "_class.txt")
    genes_status <- read.csv(pathToStatus, sep="\t")
  }, warning = function(war) {
    print("Classification file not found")
    return()
  })
  geneSt <- merge(convertionTable, genes_status, by.x="WikiGene.Name", by.y="Gene")
  return(geneSt)
}

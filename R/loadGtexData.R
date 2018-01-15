#' Coverts GTEx data from 3 files (downloaded from GTEx website) to ExpressionSet format
#'
#' This function loads files with GTEx subjects phenotypes, attributes and RNA-seq counts
#' and coverts them to one ExpressionSet objects.
#'
#'
#' @param path_to_files Path to the input folder contatining all 3 files
#' @return ExpressionSet
#' @export
load_Gtex_data <- function (path_to_files = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data")
{
  phenoFile <- paste0(path_to_files, "/GTEx_v7_Annotations_SampleAttributesDS.txt")
  pheno2File <- paste0(path_to_files, "/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
  geneFile <- paste0(path_to_files, "/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz")
  pd <- read_delim(phenoFile, delim = "\t")
  pd <- as.matrix(pd)
  rownames(pd) <- pd[, "SAMPID"]
  ids <- sapply(strsplit(pd[, "SAMPID"], "-"), function(i) paste(i[1:2], collapse = "-"))
  pd2 <-  read_delim(pheno2File, delim = "\t")
  pd2 <- as.matrix(pd2)
  rownames(pd2) <- pd2[, "SUBJID"]
  pd2 <- pd2[which(rownames(pd2) %in% unique(ids)), ]
  pd2 <- pd2[match(ids, rownames(pd2)), ]
  rownames(pd2) <- colnames(counts)
  pdfinal <- AnnotatedDataFrame(data.frame(cbind(pd, pd2)))
  cnts <- read_delim(geneFile, delim = "\t", skip = 2)
  genes <- unlist(cnts[, 1])
  geneNames <- unlist(cnts[, 2])
  counts <- cnts[, -c(1:2)]
  counts <- as.matrix(counts)
  rownames(counts) <- genes
  for (i in 1:nrow(problems(cnts))) {
    counts[problems(cnts)$row[i], problems(cnts)$col[i]] <- 1e+05
  }
  throwAway <- which(rowSums(counts) == 0)
  counts <- counts[-throwAway, ]
  genes <- sub("\\\\..*", "", rownames(counts))
  host <- "grch37.ensembl.org"
  biomart <- "ENSEMBL_MART_ENSEMBL"
  dataset <- "hsapiens_gene_ensembl"
  attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
  pdfinal <- pdfinal[match(colnames(counts), rownames(pdfinal)),]
  es <- ExpressionSet(as.matrix(counts))
  phenoData(es) <- pdfinal
  pData(es)["GTEX-YF7O-2326-101833-SM-5CVN9", "SMTS"] <- "Skin"
  pData(es)["GTEX-YEC3-1426-101806-SM-5PNXX", "SMTS"] <- "Stomach"
  es <- yarn::annotateFromBiomart(obj = es, genes = genes, host = host, biomart = biomart, dataset = dataset, attributes = attributes)
  return(es)
}

#' Subsets GTEx ExpressionSet to a given tissue
#'
#'
#'
#' @param gtex_es ExpressionSet with GTEx data
#' @param name Full tissue name to subset to
#' @return ExpressionSet
#' @export
subsetGtexData <- function (gtex_es, name)
{
  subset_tissue <- gtex_es[, which(pData(gtex_es)[, "SMTSD"] %in% name)]
  return(subset_tissue)
}


#' Saves selected tissues counts to .Rda files
#'
#'
#' @param path_to_raw_files Path to raw GTEx files
#' @param path_outFolder Path to folder to save Rda files for tissues
#' @param number_threshold Minimum number of samples per tissue to consider the tissue
#' @return
#' @export
save_Gtex_data_by_tissue <- function (path_to_raw_files = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data",
                                  path_outFolder = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/data_Rda/",
                                  samples_thr = 30)
{
  gtex_es <- load_Gtex_data()
  tissues <- names(table(pData(gtex_es)$SMTSD)[table(pData(gtex_es)$SMTSD)>=samples_thr])
  for (name in tissues) {
    subset_tissue <- subsetGtexData(gtex_es, name)
    pathSave <- paste(path_outFolder, name, sep="")
    pathSave <- paste(pathSave, ".Rda", sep="")
    save(subset_tissue, file=pathSave)
  }
}

# Run to save .Rda files to the default folder
#save_Gtex_data_by_tissue()

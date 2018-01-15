#' Plot ratio densities for all tissues
#'
#' @return ggplot2
#
plot_density_ratio_per_individual <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  for (fullName in tissues) {
    pl <- plot_density_ratio_per_individual(fullName)
    cowplot::save_plot(pl, file=paste0("results/Figures/ratio_per_individual/", fullName, "_individual_ratio_density.pdf"))
  }
}
#' Reads ratios from the file and makes a plot of density, colored by individual
#'
#' @param fullName full name of the tissue
#' @return ggplot2
#
plot_density_ratio_per_individual <- function(fullName) {
  df <- read_ratios(fullName)
  df_melted <- reshape2::melt(df)
  df_melted <- df_melted[!is.na(df_melted$value), ]
  pl <- ggplot(df_melted, aes(x = value, col = variable)) +
    geom_density(size=0.2) +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle(paste0(fullName, " #samples=", dim(df)[2]-1))
  return(pl)
}

#' Reads ratios from the file, subsamples and makes a plot
#'
#' @return
#
get_samples_tissues <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  ratios_all <- data.frame(GENE_ID = character(), ratio = numeric(), tissue = character())
  for (fullName in tissues) {
    ratios_all <- rbind(ratios_all, get_ratios_sample(fullName))
  }
  ggplot(ratios_all, aes(x = ratio, col = tissue)) +
    geom_density() +
    theme_bw()
}

#' Saves counts to files for the tissues of interest
#'
#' @return
#
save_counts_tissues_of_interest <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  for (fullName in tissues) {
    save_counts(fullName)
  }
}



#' Saves ratios to files for the tissues of interest
#'
#' @return
#
save_ratios_tissues_of_interest <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  for (fullName in tissues) {
    save_ratios(fullName)
  }
}

#' Gets ratios for a given tissue
#'
#' @param fullName full tissue name
#' @param n_samples number of rows to get
#' @return dataframe
#
get_ratios_sample <- function(fullName, n_samples = 10000) {
  df <- read_ratios(fullName)
  df_melt <- reshape2::melt(df)
  colnames(df_melt)[3] <- "ratio"
  df_melt <- df_melt[!is.na(df_melt$ratio), ]
  df_melt_sample <- dplyr::sample_n(df_melt, n_samples) %>%
    mutate(tissue = fullName) %>%
    dplyr::select(GENE_ID, ratio, tissue)
  return(df_melt_sample)
}

#' Reads and saves ratios for a given tissue
#'
#' @param fullName full tissue name
#' @param path path to save the tissues' ratios
#' @return
#
save_ratios <- function(fullName, path = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/processed_data/tissue_AI_ratios/") {
  df <- read_counts_and_calc_ratios(fullName)
  file_path <- paste0(path, fullName)
  file_path <- paste0(file_path, "_ratios.txt")
  write.table(df, file = file_path, sep = "\t", row.names = F, quote = F)
}

#' Reads ratios for a given tissue
#'
#' @param fullName full tissue name
#' @param path path with the tissues' ratios
#' @return
#
read_ratios <- function(fullName, path = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/processed_data/tissue_AI_ratios/") {
  file_path <- paste0(path, fullName)
  file_path <- paste0(file_path, "_ratios.txt")
  df <- read_delim(file_path, delim = "\t")
  return(df)
}

#' Reads and saves counts for a given tissue
#'
#' @param fullName full tissue name
#' @param path path to save the tissues' ratios
#' @return
#
save_counts <- function(fullName, path = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/processed_data/tissue_AI_counts/") {
  df <- read_counts_and_merge_samples(fullName)
  file_path <- paste0(path, fullName)
  file_path <- paste0(file_path, "_counts.txt")
  write.table(df, file = file_path, sep = "\t", row.names = F, quote = F)
}

#' Reads counts for a given tissue
#'
#' @param fullName full tissue name
#' @param path path with the tissues' ratios
#' @return
#
read_counts <- function(fullName, path = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/processed_data/tissue_AI_counts/") {
  file_path <- paste0(path, fullName)
  file_path <- paste0(file_path, "_counts.txt")
  df <- read_delim(file_path, delim = "\t")
  return(df)
}

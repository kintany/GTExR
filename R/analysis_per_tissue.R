#' Get allelic ratio's standard deviation for a tissue, individual-wise
#' @param fullName full tissue name
#' @export
#' @return dataframe
#
get_ratio_sd <- function(fullName) {
  # Read ratios
  df <- read_ratios(fullName)
  # Make the data tidy
  df_melted <- reshape2::melt(df)
  # Exclude NAs
  df_melted <- df_melted[!is.na(df_melted$value), ]
  # Calculate SDs for each individual
  df_melted_sd <- df_melted %>%
    group_by(variable) %>%
    mutate_at(.vars=vars(value), .funs = funs(sd="sd")) %>%
    select(variable, sd) %>%
    unique()
  # Convert ids
  sp <- function(x) {
    x <- gsub("ratio_", "", x)
    x <- gsub("\\.", "-", x)
  }
  df_melted_sd$id <- purrr::map_chr(df_melted_sd$variable, sp)
  df_melted_sd <- df_melted_sd[,c(3,2)]
  return(df_melted_sd)
}

#' Analyze sd distributions for all individuals for all tissues
#' @export
#' @return ggplot2
#
analyze_sd_for_tissues <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  df <- data.frame(id=character(), sd=numeric(), tissue=character())
  df_table <- data.frame(tissue=character(), fraction=character())
  for (fullName in tissues) {
    df_melted_sd <- get_ratio_sd(fullName)
    df_melted_sd$tissue <- fullName
    df <- rbind(df, df_melted_sd)
    biased <- ifelse(df_melted_sd$sd>=quantile(df_melted_sd$sd, 0.95), 1, 0)
    if (dim(table(biased))==1) tmp <- data.frame(tissue=fullName, fraction=0)
     else tmp <- data.frame(tissue=fullName, fraction=table(biased)[2] / (table(biased)[1]+table(biased)[2]))
    df_table <- rbind(df_table, tmp)
  }
  View(df_table)
  ggplot(df, aes(x=sd, color=tissue)) + geom_density()
}


#' Do variation analysis for samples with moderate sd values
#' @param fullName full tissue name
#' @export
#' @return None
variation_analysis_for_subset <- function(fullName) {
  df_melted_sd <- get_ratio_sd(fullName)
  tissues_list <- df_melted_sd$id[df_melted_sd$sd <= quantile(df_melted_sd$sd, 0.95)]
  p <- list()
  for (i in 1:6) {
    df <- process_to_counts(fullName, keepCounts = TRUE, sample_names = sample(tissues_list,2))
    colnames(df)[4:5] <- c("sample1", "sample2")
    p[[i]] <- ggplot(df, aes(x=sample1, y=sample2, color=Class)) + geom_point(size=0.5) + scale_x_log10() + scale_y_log10() + theme(legend.position = "None")
  }
  gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], ncol=3)

  df_melted_sd <- get_ratio_sd(fullName)
  tissues_list <- df_melted_sd$id[df_melted_sd$sd <= quantile(df_melted_sd$sd, 0.8)]
  for (i in 1:10) {
    df <- process_to_counts(fullName, sample_names = sample(tissues_list,10))
    df <- df[order(df$CV2),]
    top_50 <- df$Ensembl.Gene.ID[1:100]
    bottom_50 <- df$Ensembl.Gene.ID[(length(df$Ensembl.Gene.ID)-99):length(df$Ensembl.Gene.ID)]
    t <- intersect(t, top_50)
    b <- intersect(b, bottom_50)
  }
  plot_CV2_vs_mean(df, "exclude_high_sd", plotMode = "Merged", modeBoxPlot = "Sizes", n_bins = 5)
}


#' Analyze rin vs. sd for all individuals for a given tissue
#' @param fullName full tissue name
#' @export
#' @return ggplot2
#
analyze_rin <- function(fullName) {
  # Get ratios sd
  df_melted_sd <- get_ratio_sd(fullName)
  # Get RINs
  sv <- load_counts(fullName)
  rin <- data.frame(sv$SMRIN)
  rin$id <- row.names(rin)
  colnames(rin)[1] <- "RIN"
  # Merge sd and RINs
  df <- merge(df_melted_sd, rin, by.x="id", by.y="id")
  df$RIN <- as.numeric(as.character(df$RIN))
  pl <- ggplot(df, aes(x=sd, y=RIN)) +
    geom_point(size=0.5) +
    ggtitle(paste(fullName, " corr=", round(cor(df3$sd, df3$RIN),2))) +
    xlim(c(0,0.3))
  return(pl)
}

#' Produce plot with RIN vs SD(allelic_ratio) for several tissues
#' @export
#' @return ggplot2
#
get_rin_correlations_for_tissues <- function() {
  p1 <- analyze_rin("Lung")
  p2 <- analyze_rin("Liver")
  p3 <- analyze_rin("Heart - Left Ventricle")
  p4 <- analyze_rin("Spleen")
  p5 <- analyze_rin("Kidney - Cortex")
  p6 <- analyze_rin("Pancreas")
  plot_corr <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)
  cowplot::save_plot(plot_corr, file="results/Figures/corr.pdf", base_height = 10)
}

#' Outputs sample ids for individuals with high or low ratios variability
#'
#' @param fullName full tissue name
#' @param option "high" or "low" - output sampl names with high or low ratios variability
#' @export
#' @return vector with sample ids
#
get_individuals_subset <- function(fullName, option) {
  # Get ratios sd
  df_melted_sd <- get_ratio_sd(fullName)
  if (option=="high") {
    ind <- as.vector(t(df_melted_sd[df_melted_sd$sd>quantile(df_melted_sd$sd, 0.99),"id"]))
  }
  else if (option=="low") {
    ind <- as.vector(t(df_melted_sd[df_melted_sd$sd<quantile(df_melted_sd$sd, 0.01),"id"]))
  }
  else {
    print("Unknown option")
    return(NA)
  }
  return(ind)
}

#' Compare same individuals in different tissues
#' @param sample_ids ids for samples of interest
#' @export
#' @return ggplot2
#
compare_same_ind_diff_tissues <- function(sample_ids) {
  tissues <- c("Lung", "Stomach", "Small Intestine - Terminal Ileum", "Whole Blood")
  tissue_name <- purrr::map_chr(tissues, get_short_name)
  sv_lung <- process_to_counts(tissue_name[1], keepCounts = TRUE, sample_names = sample_ids)
  sv_stomach <- process_to_counts(tissue_name[2], keepCounts = TRUE, sample_names = sample_ids)
  sv_sm_int <- process_to_counts(tissue_name[3], keepCounts = TRUE, sample_names = sample_ids)
  p1 <- ggplot(sv_lung, aes(x=sv_lung[,4], y=sv_lung[,5], color=Class)) + geom_point(size=0.5) + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
  p2 <- ggplot(sv_stomach, aes(x=sv_stomach[,4], y=sv_stomach[,5], color=Class)) + geom_point(size=0.5) + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
  p3 <- ggplot(sv_sm_int, aes(x=sv_sm_int[,4], y=sv_sm_int[,5], color=Class)) + geom_point(size=0.5) + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
  gridExtra::grid.arrange(p1,p2,p3,ncol=1)
  sv_lung <- process_to_counts(tissue_name[1], keepCounts = TRUE, sample_names = c("GTEX-14BMU-0526-SM-73KW4", "GTEX-QMR6-1926-SM-32PL9"))
  ggplot(sv_lung, aes(x=sv_lung[,4], y=sv_lung[,5], color=Class)) + geom_point(size=0.5) + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
  }

#' Compares sample with high and low ratio variability
#'
#' @param fullName full tissue name
#' @export
#' @return ggplot2
#
compare_high_vs_low_samples <- function(fullName) {
  name <- get_short_name(fullName)
  # Get sample ids with low and high ratios variability
  ind_high <- get_individuals_subset(fullName, "high")
  ind_low <- get_individuals_subset(fullName, "low")
  ind <- c(ind_high, ind_low)
  # Load and subset data
  dataAll <- process_to_counts(name, keepCounts = TRUE, sample_names = ind)
  # Pairwise comparison: 2 samples with low ratios variability
  d1 <- process_to_counts(name, keepCounts = TRUE, sample_names = ind_low[3:4])
  # And two samples with high
  d2 <- process_to_counts(name, keepCounts = TRUE, sample_names = ind_high[1:2])
  # Merge
  d12 <- merge(d1[,c(1,3,10)], d2[,c(1,10)], by.x="Ensembl.Gene.ID", by.y="Ensembl.Gene.ID")
  colnames(d12) <- c("Ensembl.Gene.ID","Class","CV2.samples12","CV2.samples34" )
  ggplot(d12, aes(x=CV2.samples12, y=CV2.samples34, color=Class)) +
    geom_point(size=0.3) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
  d12_melted <- reshape2::melt(d12)
  d12_melted$id_u <- paste0(d12_melted$Class, "_", d12_melted$variable)
  ggplot(d12_melted, aes(y=value, x=id_u, fill=id_u)) +
    geom_boxplot() +
    scale_y_log10()
}


#' Plot ratio densities for all tissues
#' @export
#' @return ggplot2
#
plot_density_ratio_per_individual_all_tissues <- function() {
  tissues <- c(
    "Kidney - Cortex", "Heart - Left Ventricle",
    "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen"
  )
  for (fullName in tissues) {
    pl <- plot_density_ratio_per_individual(fullName, byStatus = TRUE)
    cowplot::save_plot(pl, file=paste0("results/Figures/ratio_per_individual/", fullName, "_individual_ratio_density_MAE.pdf"))
  }

  pltList <- list()
  y <- 0
  for (fullName in tissues) {
    pl <- plot_density_ratio_per_individual(fullName, byStatus = FALSE)
    pltList[[y+1]] <- pl
    y <- y+1
  }
  pl <- gridExtra::marrangeGrob(pltList, ncol=2, nrow=4)
  cowplot::save_plot(pl, file="results/Figures/ratio_per_individual/all_tissues.pdf", base_height = 12)
}

#' Reads ratios from the file and makes a plot of density, colored by individual
#'
#' @param fullName full name of the tissue
#' @param byStatus if TRUE, make plots for MAE/BAE status
#' @export
#' @return ggplot2
#
plot_density_ratio_per_individual <- function(fullName, byStatus = FALSE) {
  df <- read_ratios(fullName)
  df_melted <- reshape2::melt(df)
  df_melted <- df_melted[!is.na(df_melted$value), ]
  if (byStatus) {
    short_name <- get_short_name(fullName)
    geneSt <- get_classification(short_name)
    df_melted_class <- merge(df_melted, geneSt, by.x="GENE_ID", by.y="Ensembl.Gene.ID")
    pl1 <- ggplot(df_melted_class, aes(x = value, col = Class, linetype = variable)) +
      geom_density(size=0.2) +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle(paste0(fullName, " #samples=", dim(df)[2]-1)) +
      scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
    pl2 <- ggplot(df_melted_class, aes(x = value, col = Class)) +
      geom_density(size=0.2) +
      theme_bw() +
      ggtitle(paste0(fullName, " #samples=", dim(df)[2]-1)) +
      scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group")
    pl <- gridExtra::grid.arrange(pl1, pl2, nrow = 2)
  }
  else {
    pl <- ggplot(df_melted, aes(x = value, col = variable)) +
      geom_density(size=0.2) +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle(paste0(fullName, " #samples=", dim(df)[2]-1))
  }
  return(pl)
}

#' Reads counts, apply thresholds, plot distributions
#'
#' @param fullName full tissue name
#' @export
#' @return ggplot2
#
analyze_sd_by_threshold <- function(fullName) {
  # Read counts
  counts <- read_counts(fullName)
  # Calculate ratios
  counts <- counts %>%
    mutate(ratio = ref_count / (ref_count + alt_count))
  counts_thresholds <- c(8, 16, 32, 64)
  y <- 1
  pltList <- list()
  for (thr in counts_thresholds) {
    # Subset counts
    counts_thr <- counts %>% filter(ref_count + alt_count >= thr)
    pltList[[y]] <- ggplot(counts_thr, aes(x = ratio, col = sample_id)) +
      geom_density(size=0.2) +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle(paste0("thresold=", thr))
    y <- y + 1
  }
  gridExtra::grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], ncol=2)
  #p <- gridExtra::grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], ncol=2)
  #cowplot::save_plot(p, file="results/Figures/Lung_sd_by_threshold.pdf", base_height = 8)
}

#' Reads ratios from the file and analyze standard deviations
#'
#' @param fullName full tissue name
#' @export
#' @return ggplot2
#
analyze_sd_MAE_BAE <- function(fullName) {
  # Read counts
  counts <- read_counts(fullName)
  # Calculate ratios
  counts_1 <- counts %>%
    mutate(ratio = ref_count / (ref_count + alt_count))
  # Merge with classification
  short_name <- get_short_name(fullName)
  geneSt <- get_classification(short_name)
  df <- merge(counts_1, geneSt, by.x="GENE_ID", by.y="Ensembl.Gene.ID")
  # Do the analysis for several expression bins
  pltList <- list()
  y <- 1
  counts_thresholds = c(10, 50, 100, 150, 200, 10000)
  for (thr in counts_thresholds) {
    # Subset counts
    df_thr <- df %>% filter(ref_count + alt_count >= counts_thresholds[y] &  ref_count + alt_count < counts_thresholds[y+1])
    # Calculate SD
    df1 <- df_thr %>%
      group_by(Class, sample_id) %>%
      mutate_at(.vars=vars(ratio), .funs = funs(sd="sd")) %>%
      select(sample_id, Class, sd) %>%
      unique()
    # Make a plot
    a <- df_thr %>% select(6,7) %>% unique()
    pltList[[y]] <- ggplot(df1, aes(x=Class, y=sd, fill=Class)) +
      geom_boxplot() +
      theme(legend.position = "None") +
      ylim(c(0,0.3)) +
      ggtitle(paste0("# MAE=", table(a$Class)[2], " # BAE=", table(a$Class)[1]))
    y <- y + 1
  }
  p <- gridExtra::grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], pltList[[6]], ncol=3)
  cowplot::save_plot(p, file="results/Figures/sd_MAE_BAE.pdf", base_height = 8)
}

#' Reads ratios from the file and calculate gini indexes
#'
#' @param fullName full tissue name
#' @export
#' @return ggplot2
#
analyze_imbalance_vs_counts <- function(fullName) {
  ratios <- read_ratios(fullName)
  name <- get_short_name(fullName)
  df <- process_to_counts(name, keepCounts = TRUE)
  df_cor <- data.frame(gene_id=character(), cor=numeric(), c=numeric())
  ratios <- data.frame(ratios)
  cor_an <-function(x, gene = FALSE) {
    gene_id <- as.character(x[1])
    counts_string <- data.frame(t(df[df$Ensembl.Gene.ID==gene_id, 4:(dim(df)[2]-6)]))
    if (dim(counts_string)[2]>0) {
      counts_string$gene_id <- row.names(counts_string)
      colnames(counts_string)[1] <- "counts"
      x <- x[2:length(x)]
      x <- data.frame(t(x))
      colnames(x) <- "ratio"
      x$gene_id <- row.names(x)
      x <- x[!is.na(x$ratio),]
      sp <- function(x) {
        x <- gsub("ratio_", "", x)
      }
      x$gene_id <- purrr::map_chr(x$gene_id, sp)
      x$ratio_abs <- abs(x$ratio-0.5)
      merge_x <- merge(x, counts_string, by.x="gene_id", by.y="gene_id")
      pl <- ggplot(merge_x, aes(x=ratio_abs, y=counts)) +
        geom_point() +
        ggtitle(gene_id) +
        geom_smooth(method='lm')
      ggplot(merge_x, aes(x=ratio, y=counts)) +
        geom_point() +
        ggtitle(gene_id)
      tmp_cor <- data.frame(gene_id=gene_id, cor=cor(merge_x$ratio_abs, merge_x$counts), c=dim(merge_x)[1])
      if (gene) return(pl)
      else return(tmp_cor)
    }
    else {
      return(NA)
    }
  }
  y <- 1
  tmp_cor_gene <- list()
  for (i in 1:dim(ratios)[1]) {
    tmp_cor <- cor_an(ratios[i,])
    if (!is.na(tmp_cor)) {
      df_cor <- rbind(df_cor, tmp_cor)
      if (tmp_cor$cor<= -0.6 & tmp_cor$c > 20) {
        tmp_cor_gene[[y]] <- cor_an(ratios[i,], gene = T)
        y <- y+1
      }
    }
  }
  p <- gridExtra::grid.arrange(tmp_cor_gene[[1]], tmp_cor_gene[[2]], tmp_cor_gene[[3]],
                               tmp_cor_gene[[4]], tmp_cor_gene[[5]], tmp_cor_gene[[6]],
                               tmp_cor_gene[[7]], tmp_cor_gene[[8]], ncol=4)
  cowplot::save_plot(p, file="results/Figures/counts_vs_ratio_exemples_Kidney.pdf", base_height = 8, base_aspect_ratio = 2)
}
#' Run ratios analysis for NPCs
#' @export
#' @return dataframe
#
make_plot_for_NPC <- function() {
  G_data <- read.delim("~/Dropbox (Partners HealthCare)/variation_project/Gendrel_data/G_data.txt")
  G_data_melt <- reshape2::melt(G_data)
  ggplot(G_data_melt, aes(x=value, color=res)) + geom_density(adjust = 3)
  G_data_melt_2 <- G_data_melt
  G_data_melt_2$variable <- "mix"
  G_data_melt_2 <- G_data_melt_2 %>%
    group_by(ensembl_gene_id) %>%
    mutate_at(.vars=vars(value), .funs = funs(mean="mean")) %>%
    select(ensembl_gene_id, variable, mean) %>%
    unique()
  G_data_melt_2 <- data.frame(G_data_melt_2)
  G_data_melt <- G_data_melt[,c(1,3,4)]
  colnames(G_data_melt_2) <- colnames(G_data_melt)
  G <- rbind(G_data_melt, G_data_melt_2)
  ggplot(G, aes(x=value, color=variable)) + geom_density(adjust = 3) + scale_color_manual(values=c(rep("blue", 8), "red"))
  G <- G[!is.na(G$value),]
  G_sd <- G %>%
    group_by(variable) %>%
    mutate_at(.vars=vars(value), .funs = funs(sd="sd")) %>%
    select(variable, sd) %>%
    unique()
  G_sd$ggini <- c(apply(G_data[,c(3:10)],2,edgeR::gini), edgeR::gini(G$value[G$variable=="mix"]))
  retirn(G_sd)
}

#' Run ratios analysis for NPCs
#' @export
#' @param sample_id sample id
#' @return dataframe
#
list_of_samples_with_ratios <- function(sample_id) {
  # Get list of available tissues for a given individual
  path_to_files <-  "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data"
  phenoFile <- paste0(path_to_files, "/GTEx_v7_Annotations_SampleAttributesDS.txt")
  pd <- read_delim(phenoFile, delim = "\t")
  pd_sub <- pd[startsWith(pd$SAMPID, sample_id),]
  pd_sub <- pd_sub %>%
    filter(SMRIN >= 6) %>%
    dplyr::select(SAMPID, SMTSD)
  # Get tissue codes
  pd_sub$code <- purrr::map_chr(pd_sub$SMTSD, get_tissue_code)
  # Read file with ASE data for a given individual
  f_name <- paste0("/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/data/raw_data/phe000024.v1.GTEx_ASE_SNPs.expression-matrixfmt-ase.c1/", sample_id)
  f_name <- paste0(f_name, ".ase_table.tsv.gz")
  df <- readr::read_delim(f_name, "\t", escape_double = FALSE, trim_ws = TRUE)
  # Subset and calculate ratios SNP-wise
  df <- df %>% filter(TISSUE_ID %in% pd_sub$code) %>% filter(LOW_MAPABILITY ==0 & MAPPING_BIAS_SIM == 0 & GENOTYPE_WARNING ==0)
  df <- df[,c(21,1:6,8,9,10,19,20)]
  df <- df[!is.na(df$GENE_ID),]
  df$sw <- ifelse(df[,11]== "GT;1|0", 1, 0)
  df$REF_COUNT_UPD <- df$REF_COUNT
  df$ALT_COUNT_UPD <- df$ALT_COUNT
  df[df$sw==1,"REF_COUNT_UPD"] <- df[df$sw==1,"ALT_COUNT"]
  df[df$sw==1,"ALT_COUNT_UPD"] <- df[df$sw==1,"REF_COUNT"]
  #df$RATIO <- df$REF_COUNT_UPD / (df$REF_COUNT_UPD + df$ALT_COUNT_UPD)
  # Aggregate by tissue
  res <- df %>%
    group_by(SAMPLE_ID, GENE_ID) %>%
    mutate_at(.vars=vars(REF_COUNT_UPD, ALT_COUNT_UPD), .funs = funs(sum="sum")) %>%
    select(GENE_ID, TISSUE_ID, REF_COUNT_UPD_sum, ALT_COUNT_UPD_sum) %>%
    unique() %>%
    filter(REF_COUNT_UPD_sum + ALT_COUNT_UPD_sum >= 10) %>%
    mutate(RATIO = REF_COUNT_UPD_sum / (REF_COUNT_UPD_sum + ALT_COUNT_UPD_sum)) %>%
    group_by(TISSUE_ID) %>%
    mutate_at(.vars=vars(RATIO), .funs = funs(sd="sd")) %>%
    select(SAMPLE_ID, TISSUE_ID, sd) %>%
    unique()
  return(res)
}


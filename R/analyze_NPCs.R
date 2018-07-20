#' Runs analysis for NPCs
#'
#' @export
#' @return ggplot
#
analyze_NPC_var <- function() {
  Table_S1_mod <- read.delim("~/Dropbox (Partners HealthCare)/variation_project/Gendrel_data/data/Table_S1_mod_minus_chrX_parsed_mod.txt")
  Table_S1_mod <- Table_S1_mod %>% group_by(Gene) %>% slice(which.max(GEI.71b_counts))
  NPCs_counts <- Table_S1_mod[,c(1:18)]
  NPCs_counts <- NPCs_counts %>% filter(NPC_Catalog %in% c("Biallelic", "Monoallelic")) %>% as.data.frame
  dataNorm <- get_dataNorm_from_matrix(NPCs_counts[,c(1:10)])
  orPlot <- plot_dataNorm(dataNorm)
  save_plot(orPlot, file="results/NPCs/8clones_plot.pdf", base_width = 14)

  #All possible combination of pairs
  combin <- combn(1:8, 2)
  pltList <- list()
  for (j in 1:6) {
    NPCs_counts_merged <- selectPairs(NPCs_counts, combin[,sample(1:28, 8)])
    dataNorm <- get_dataNorm_from_matrix(NPCs_counts_merged)
    pltList[[j]] <- plot_dataNorm(dataNorm)
  }
  pair_plot <- grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], pltList[[6]], ncol=3)
  save_plot(pair_plot, file="results/NPCs/pairs.pdf", base_width = 14)

  #All possible combination of trios
  combin <- combn(1:8, 3)
  pltList <- list()
  for (j in 1:6) {
    NPCs_counts_merged <- selectTrios(NPCs_counts, combin[,sample(1:56, 8)])
    dataNorm <- get_dataNorm_from_matrix(NPCs_counts_merged)
    pltList[[j]] <- plot_dataNorm(dataNorm)
  }
  pair_plot <- grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], pltList[[6]], ncol=3)
  save_plot(pair_plot, file="results/NPCs/trios.pdf", base_width = 14)

  #All possible combination of 4
  combin <- combn(1:8, 4)
  pltList <- list()
  for (j in 1:6) {
    NPCs_counts_merged <- select4(NPCs_counts, combin[,sample(1:70, 8)])
    dataNorm <- get_dataNorm_from_matrix(NPCs_counts_merged)
    pltList[[j]] <- plot_dataNorm(dataNorm)
  }
  pair_plot <- grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], pltList[[6]], ncol=3)
  save_plot(pair_plot, file="results/NPCs/comb4.pdf", base_width = 14)

  #All possible combination of 5
  combin <- combn(1:8, 5)
  pltList <- list()
  for (j in 1:6) {
    NPCs_counts_merged <- select5(NPCs_counts, combin[,sample(1:56, 8)])
    dataNorm <- get_dataNorm_from_matrix(NPCs_counts_merged)
    pltList[[j]] <- plot_dataNorm(dataNorm)
  }
  pair_plot <- grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], pltList[[6]], ncol=3)
  save_plot(pair_plot, file="results/NPCs/comb5.pdf", base_width = 14)
}

selectPairs <- function(NPCs_counts, x) {
  NPCs_counts %>%
    mutate(m1 = (.[[x[1,1]+2]] + .[[x[2,1]+2]])/2,
           m2 = (.[[x[1,2]+2]] + .[[x[2,2]+2]])/2,
           m3 = (.[[x[1,3]+2]] + .[[x[2,3]+2]])/2,
           m4 = (.[[x[1,4]+2]] + .[[x[2,4]+2]])/2,
           m5 = (.[[x[1,5]+2]] + .[[x[2,5]+2]])/2,
           m6 = (.[[x[1,6]+2]] + .[[x[2,6]+2]])/2,
           m7 = (.[[x[1,7]+2]] + .[[x[2,7]+2]])/2,
           m8 = (.[[x[1,8]+2]] + .[[x[2,8]+2]])/2) %>%
    select(Gene, NPC_Catalog, m1, m2, m3, m4, m5, m6, m7, m8)
}

selectTrios <- function(NPCs_counts, x) {
  NPCs_counts %>%
    mutate(m1 = (.[[x[1,1]+2]] + .[[x[2,1]+2]] + .[[x[3,1]+2]])/3,
           m2 = (.[[x[1,2]+2]] + .[[x[2,2]+2]] + .[[x[3,2]+2]])/3,
           m3 = (.[[x[1,3]+2]] + .[[x[2,3]+2]] + .[[x[3,3]+2]])/3,
           m4 = (.[[x[1,4]+2]] + .[[x[2,4]+2]] + .[[x[3,4]+2]])/3,
           m5 = (.[[x[1,5]+2]] + .[[x[2,5]+2]] + .[[x[3,5]+2]])/3,
           m6 = (.[[x[1,6]+2]] + .[[x[2,6]+2]] + .[[x[3,6]+2]])/3,
           m7 = (.[[x[1,7]+2]] + .[[x[2,7]+2]] + .[[x[3,7]+2]])/3,
           m8 = (.[[x[1,8]+2]] + .[[x[2,8]+2]] + .[[x[3,8]+2]])/3) %>%
    select(Gene, NPC_Catalog, m1, m2, m3, m4, m5, m6, m7, m8)
}

select4 <- function(NPCs_counts, x) {
  NPCs_counts %>%
    mutate(m1 = (.[[x[1,1]+2]] + .[[x[2,1]+2]] + .[[x[3,1]+2]] + .[[x[4,1]+2]])/4,
           m2 = (.[[x[1,2]+2]] + .[[x[2,2]+2]] + .[[x[3,2]+2]] + .[[x[4,2]+2]])/4,
           m3 = (.[[x[1,3]+2]] + .[[x[2,3]+2]] + .[[x[3,3]+2]] + .[[x[4,3]+2]])/4,
           m4 = (.[[x[1,4]+2]] + .[[x[2,4]+2]] + .[[x[3,4]+2]] + .[[x[4,4]+2]])/4,
           m5 = (.[[x[1,5]+2]] + .[[x[2,5]+2]] + .[[x[3,5]+2]] + .[[x[4,5]+2]])/4,
           m6 = (.[[x[1,6]+2]] + .[[x[2,6]+2]] + .[[x[3,6]+2]] + .[[x[4,6]+2]])/4,
           m7 = (.[[x[1,7]+2]] + .[[x[2,7]+2]] + .[[x[3,7]+2]] + .[[x[4,7]+2]])/4,
           m8 = (.[[x[1,8]+2]] + .[[x[2,8]+2]] + .[[x[3,8]+2]] + .[[x[4,8]+2]])/4) %>%
    select(Gene, NPC_Catalog, m1, m2, m3, m4, m5, m6, m7, m8)
}

select5 <- function(NPCs_counts, x) {
  NPCs_counts %>%
    mutate(m1 = (.[[x[1,1]+2]] + .[[x[2,1]+2]] + .[[x[3,1]+2]] + .[[x[4,1]+2]] + .[[x[4,1]+2]])/5,
           m2 = (.[[x[1,2]+2]] + .[[x[2,2]+2]] + .[[x[3,2]+2]] + .[[x[4,2]+2]] + .[[x[4,2]+2]])/5,
           m3 = (.[[x[1,3]+2]] + .[[x[2,3]+2]] + .[[x[3,3]+2]] + .[[x[4,3]+2]] + .[[x[4,3]+2]])/5,
           m4 = (.[[x[1,4]+2]] + .[[x[2,4]+2]] + .[[x[3,4]+2]] + .[[x[4,4]+2]] + .[[x[4,4]+2]])/5,
           m5 = (.[[x[1,5]+2]] + .[[x[2,5]+2]] + .[[x[3,5]+2]] + .[[x[4,5]+2]] + .[[x[4,5]+2]])/5,
           m6 = (.[[x[1,6]+2]] + .[[x[2,6]+2]] + .[[x[3,6]+2]] + .[[x[4,6]+2]] + .[[x[4,6]+2]])/5,
           m7 = (.[[x[1,7]+2]] + .[[x[2,7]+2]] + .[[x[3,7]+2]] + .[[x[4,7]+2]] + .[[x[4,7]+2]])/5,
           m8 = (.[[x[1,8]+2]] + .[[x[2,8]+2]] + .[[x[3,8]+2]] + .[[x[4,8]+2]] + .[[x[4,8]+2]])/5) %>%
    select(Gene, NPC_Catalog, m1, m2, m3, m4, m5, m6, m7, m8)
}

#' Plots CV2 vs mean
#' @param input data frame with all data
#' @export
#' @return ggplot
#
plot_dataNorm <- function(df) {
  pl <- ggplot(df, aes(x=rowMean, y=CV2, color=factor(NPC_Catalog))) +
    geom_point(aes(alpha = 0.5), size=0.5) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Mean expression, normalized counts") +
    ylab("CV^2") +
    stat_smooth(aes(group=factor(NPC_Catalog)), col=NA, method = "loess", size=1, se=TRUE) +
    stat_smooth(aes(color=factor(NPC_Catalog)), method = "loess", size=1, se=F) +
    theme(legend.justification=c(1,1), legend.position=c(0.99,0.99),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group") +
    scale_alpha(guide = 'none')
  return(pl)
}

#' Calcultes all statistics for counts
#' @param input matrix with counts
#' @export
#' @return dataframe
#
get_dataNorm_from_matrix <- function(df) {
  sv_counts <- as.matrix(df[3:dim(df)[2]])
  row.names(sv_counts) <- df$Gene
  sv_counts <- sv_counts[rowSums(sv_counts)>=10, ]
  sv_counts <- sv_counts[rowSums(sv_counts>0)>(dim(sv_counts)[2]/2), ]
  lib.size <- DESeq2::estimateSizeFactorsForMatrix(sv_counts)
  nmat <- as.matrix(t(t(sv_counts)/lib.size))
  nmat <- log(nmat)
  dataNorm <- data.frame(nmat)
  dataNorm$rowSd <- matrixStats::rowSds(nmat)
  dataNorm$rowMean <- rowMeans(nmat)
  dataNorm$rowVar <- matrixStats::rowVars(nmat)
  dataNorm$CV <- matrixStats::rowSds(nmat) / rowMeans(nmat)
  dataNorm$CV2 <- dataNorm$CV * dataNorm$CV
  dataNorm$gene <- row.names(dataNorm)
  dataNorm <- merge(df, dataNorm[,c("rowSd", "rowMean", "rowVar", "CV", "CV2", "gene")], by.x="Gene", by.y="gene")
}

#' Runs analysis for NPCs: pull clones, look at SDs
#'
#' @export
#' @return ggplot
#
analyze_NPC_pulled_clones <- function() {
  Table_S1_mod <- read.delim("~/Dropbox (Partners HealthCare)/variation_project/Gendrel_data/data/Table_S1_mod_minus_chrX_parsed_mod.txt")
  Table_S1_mod <- Table_S1_mod %>% group_by(Gene) %>% slice(which.max(GEI.71b_counts))
  NPCs_counts <- Table_S1_mod[,c(1:18)]
  NPCs_counts <- NPCs_counts %>% filter(NPC_Catalog %in% c("Biallelic", "Monoallelic")) %>% as.data.frame
  # plot sds for combinations pulled of clones
  simulate_ratios(NPCs_counts)
}


#' Makes boxplots for SDs
#'
#' @param df dataframe with counts and ratios
#' @export
#' @return ggplot
#
simulate_ratios <- function(df) {
  sd_all <- data.frame(matrix(rep(NA,560), ncol=8, nrow=70))
  colnames(sd_all) <- paste0("set_of_", 1:8)
  for (y in 1:8) {
    sd_all[, y] = get_sd_ratios(df, y)
  }
  sd_all <- sd_all %>% gather()
  ggplot(sd_all, aes(y=value, x=key)) +
    geom_boxplot() +
    xlab("combinations") +
    ylab("sd")
}

#' Calcultates SDs for allelic ratios
#'
#' @param df dataframe with counts and ratios,
#' @param y number of combinations
#' @export
#' @return ggplot
#
get_sd_ratios <- function(df, y) {
  combin <- combn(1:8, y)
  sds <- rep(NA, 70)
  for (j in 1:dim(combin)[2]) {
    a <- subset_and_get_ratio(df[,3:18], combin[,j])
    sds[j] <- sd(a$ref_sum, na.rm = T)
  }
  return(sds)
}

#' Pulls data and gets allelic ratios
#'
#' @param df dataframe with counts and ratios,
#' @param x vector with columns to pull
#' @export
#' @return ggplot
#
subset_and_get_ratio <- function(df, x) {
  ratio_column_ids <- x + 8
  df_sub <- df %>% select(x, ratio_column_ids)
  for (i in 1:length(x)) {
    df_sub[,paste0("ref",i)] <- as.numeric(as.character(df_sub[,i])) * as.numeric(as.character(df_sub[,(i+length(x))]))
    df_sub[,paste0("alt",i)] <- as.numeric(as.character(df_sub[,i])) * (1-as.numeric(as.character(df_sub[,(i+length(x))])))
  }
  ref_sum <- df_sub %>% select(starts_with("ref")) %>% mutate(ref_sum = rowSums(.)) %>% select(ref_sum)
  alt_sum <- df_sub %>% select(starts_with("alt")) %>% mutate(alt_sum = rowSums(.)) %>% select(alt_sum)
  ratio <- ref_sum / (ref_sum + alt_sum)
}

#' Calculates CVs for "observed" BAE and all genes from counts and ratios
#'
#' @param a vector, first 8 entries are counts, second 8 entries are ratios
#' @export
#' @return vector with CVs and mean
#
measure <- function(a) {
  counts <- a[1:8]
  x <- a[9:16]
  MAE_more <- vector('numeric')
  MAE_less <- vector('numeric')
  BAE_x <- vector('numeric')
  i = 1
  for (x1 in x) {
    if (!is.na(x1)) {
      if (x1<0.15) {
        MAE_less <- c(MAE_less, counts[i])
      }
      else if (x1>0.85) {
        MAE_more <- c(MAE_more, counts[i])
      }
      else if ((x1>0.33)&(x1<0.66)){
        BAE_x <- c(BAE_x, counts[i])
      }
    }
    i = i + 1
  }
  cv_true_BAE <- NA
  cv_MAE <- NA
  cv_BAE <- NA
  cv_all <- NA
  cv_all_arr <- rep(0,5)
  mean_MAE <- NA
  mean_BAE <- NA
  mean_all <- NA
  mean_all_arr <- rep(0,5)
  MAE_less <- unlist(MAE_less)
  MAE_more <- unlist(MAE_more)
  BAE_x <- unlist(BAE_x)
  if (length(BAE_x)>2) {
    if (length(MAE_less)>2) {
      cv_MAE <- sd(MAE_less) / mean(MAE_less)
      cv_BAE <- sd(BAE_x) / mean(BAE_x)
      mean_MAE <- mean(MAE_less)
      mean_BAE <- mean(BAE_x)
      #print(paste("case: BAE list ",BAE_x, "BAE_sd=",sd(BAE_x), "BAE_mean=",mean(BAE_x)))
    }
    if (length(MAE_more)>2) {
      cv_MAE <- sd(MAE_more) / mean(MAE_more)
      cv_BAE <- sd(BAE_x) / mean(BAE_x)
      mean_MAE <- mean(MAE_more)
      mean_BAE <- mean(BAE_x)
    }
  }
  for (i in 1:10) {
    counts_s <- sample(counts, length(BAE_x))
    cv_all_arr[i] <- sd(as.numeric(counts_s)) / mean(as.numeric(counts_s))
    mean_all_arr[i] <- mean(as.numeric(counts_s))
  }
  cv_all <- mean(cv_all_arr)
  mean_all <- mean(mean_all_arr)
  return(c(cv_BAE,cv_all,mean_all))
}

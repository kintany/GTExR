#' Generate variatility plots for a given tissue
#'
#' @param name tissue name
#' @param mode "norm_counts" / "log" meaning using normalized by library size counts or using log(normalized_counts)
#' @return  None
#' @export
#'
generatePlots <- function(name, mode) {
  df <- process_to_counts(name, mode="log")
  plot_CV2_vs_mean(df, name, plotMode = "Merged", modeBoxPlot = "Sizes", n_bins = 5)
}

#' Convert normalized counts to z-scores and plot heatmaps
#'
#' @param name Tissue name
#' @param sortedBy CV2 or status
#' @return  None
#' @export
#'
zscores_heatmap <- function(name, sortedBy = "CV2") {
  # Read normalized counts for a given tissue
  df <- process_to_counts(name, keepCounts = T, mode="norm_counts")
  df_mod <- df
  # Function that calculates scores from counts
  getDelta <- function(x) {
    power_value <- 2
    gene_values <- log(x)
    mean_gene <- mean(gene_values)
    sd_gene <- sd(gene_values)
    gene_values_out <- (gene_values - mean_gene) * power_value
    gene_values_out <- 1 / (1 + exp(-1 * gene_values_out))
    return(gene_values_out)
  }
  # Next we apply this function to convert normalized counts to scores
  df_mod[,4:(dim(df_mod)[2]-6)] <- t(apply(df_mod[,4:(dim(df_mod)[2]-6)], 1, getDelta))
  # Create gradient palette, so that blue is low values, red is high values
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
  #my_palette <- colorRampPalette(c("white", "blue"))(n = 1000)
  # Get some test data
  test <- df_mod[1:1000,]
  if (sortedBy == "CV2") {
    test <- test[order(test$CV2),]
    }
  if (sortedBy == "status") {
    test <- test[order(test$Class),]
  }
  # Or, we can color them by MAE
  cols <- rep('#FC6621', nrow(test))
  cols[test$Class == "MAE"] <- '#105CFB'
  rownames(test) <- test$Ensembl.Gene.ID
  # Save plot
  #pdf("results/Figures/heat_1000.pdf",width=40,height=60)
  gplots::heatmap.2(as.matrix(test[,4:100]),
            Rowv=FALSE,
            Colv=FALSE,
            col=my_palette,
            colRow = cols,
            density.info="none",
            dendrogram="none",
            trace="none",
            key=FALSE,
            lwid=c(0.1,4),
            lhei=c(0.1,4),
            margins = c(20, 8)
  )
  #dev.off()
}


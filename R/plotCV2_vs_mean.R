#' Plots CV2 versus mean
#'
#'
#' @param df dataframe with statistics
#' @param name tissue name
#' @param path_fig path to save figures
#' @param n_bins number of bins for mean expression
#' @param plotMode "Both" or "Boxplots" or "Scatterplot" or "Merged" (merge boxplot and sctterplot into one plot)
#' @param modeBoxPlot "None", "Ranges" (bin ranges will be added to boxplots), "Sizes" (sample sizes will be added to the plot) or "Pvals" (pvalues will be added to boxplots)
#' @return ggplot
#' @export
#'

plot_CV2_vs_mean <- function(df,
                            name,
                            path_fig = "/Users/svetlana/Dropbox (Partners HealthCare)/variation_project/GTEx/GTExR/results/Figures/CV2_mean/",
                            n_bins = 5,
                            plotMode = "Both",
                            modeBoxPlot = "None") {
  # Create names for the plots
  fig1Name <- paste0(path_fig, name)
  fig1Name <- paste0(fig1Name, "_fig1.pdf")
  fig2Name <- paste0(path_fig, name)
  fig2Name <- paste0(fig2Name, "_fig2.pdf")
  if (plotMode == "Merged") {
    fig3Name <- paste0(path_fig, name)
    fig3Name <- paste0(fig3Name, "_fig3.pdf")
  }
  #### Sample for 5 interquantile ranges
  qq <- rep(0,n_bins)
  # Get quantile thresholds
  for (y in 0:(n_bins-1)) {
    qq[y+1] <- quantile(df$rowMean, y/n_bins)
  }
  qq[n_bins+1] <- max(df$rowMean)
  if (plotMode != "Scatterplot") {
    # First plot: boxplots
    plot1 <- plot_one(df, n_bins, qq, modeBoxPlot)
    if (!(is.na(plot1))) {
      cowplot::save_plot(fig1Name, plot1, base_aspect_ratio = 3)
    }
    else {
      print("Some bins don't have enough genes, try decreasing the number of bins")
    }
  }
  if (plotMode != "Boxplots") {
    # Second plot: scatterplot
    plot2 <- plot_two(df, n_bins, qq)
    cowplot::save_plot(fig2Name, plot2, base_height = 8)
  }
  if (plotMode == "Merged") {
    if (!(is.na(plot1))) {
      gr <- gridExtra::grid.arrange(plot2, plot1)
      cowplot::save_plot(fig3Name, gr, base_height = 10)
    }
    else {
      print("Some bins don't have enough genes, try decreasing the number of bins")
      cowplot::save_plot(fig3Name, plot2, base_height = 8)
    }
  }
}

#' Plots boxplots for ranges
#'
#'
#' @param df dataframe with statistics
#' @param n_bins number of bins for mean expression
#' @param qq quantiles
#' @param modeBoxPlot "None", "Ranges" (bin ranges will be added to boxplots), "Sizes" (sample sizes will be added to the plot) or "Pvals" (pvalues will be added to boxplots)
#' @return ggplot
#'
plot_one <- function(df, n_bins, qq, modeBoxPlot = "None") {
  pltList <- list()
  pval <- rep(0,n_bins)
  range <- rep(0,n_bins)
  for (y in 0:(n_bins-1)) {
    # get the data for the bin
    from <- df[(df$rowMean > (qq[y+1]))&(df$rowMean <= (qq[y+2])),]
    # ranges of values within bins
    range[y+1] <- paste(round(qq[y+1],2), round(qq[y+2],2), sep="-")
    # number of genes per class
    len1 <- length(from$CV2[from$Class=="MAE"])
    len2 <- length(from$CV2[from$Class=="BAE"])
    len <- min(len1, len2)
    # sample data: 300 sample (or less if one of the classes had less entries)
    if (len < 10) {
      return(NA)
    }
    else if (len < 300) {
      n1 <- sample(from$CV2[from$Class=="MAE"], len, replace = F)
      n2 <- sample(from$CV2[from$Class=="BAE"], len, replace = F)
    }
    else {
      n1 <- sample(from$CV2[from$Class=="MAE"], 300, replace = F)
      n2 <- sample(from$CV2[from$Class=="BAE"], 300, replace = F)
    }
    n1_n <- cbind(n1, "MAE", "#FC6621")
    n2_n <- cbind(n2, "BAE", "#105CFB")
    # create dataframe with subset data
    dat_plot <- data.frame(rbind(n1_n, n2_n), stringsAsFactors = F)
    colnames(dat_plot) <- c("CV2", "Class", "Color")
    # Create boxplot for the data
    plotBox <- ggplot(dat_plot, aes(y=as.numeric(CV2), x=factor(Class), fill=factor(Class), color=factor(Class))) +
      geom_boxplot() +
      scale_fill_manual(values=c("#FC6621", "#105CFB")) +
      scale_color_manual(values=c("#FC6621", "#105CFB")) +
      xlab((y+1)) +
      ylab("") +
      scale_y_log10(limits = c(0.1, 3)) +
      theme_bw() +
      stat_summary(geom = "crossbar", width = 0.6, fatten=0, color="white",
                   fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})
    if (y!=n_bins) {
      plotBox <- plotBox +
        theme(legend.position = "none")
    }
    # Calculate p-values
    pval[y+1] <- signif(wilcox.test(n1,n2, alternative = "greater")$p.value, 4)
    # Add labels if needed
    if (modeBoxPlot == "Ranges")
      plotBox <- plotBox + annotate("text", x = 1.5, y = 2.5, label = range[y+1])
    if (modeBoxPlot == "Sizes")
      plotBox <- plotBox + annotate("text", x = 1.5, y = 2.5, label = len)
    if (modeBoxPlot == "Pvals")
      plotBox <- plotBox + annotate("text", x = 1.5, y = 2.5, label = paste0("p=", pval[y+1]))
    pltList[[y+1]] <- plotBox
  }
  #plot1 <- gridExtra::marrangeGrob(pltList, ncol=n_bins, nrow=1)
  plot1 <- gridExtra::grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], pltList[[4]], pltList[[5]], ncol=5)
  return(plot1)
}
#' Plots scatterplots
#'
#'
#' @param df dataframe with statistics
#' @param n_bins number of bins for mean expression
#' @param qq quantiles
#' @return ggplot
#'
plot_two <- function(df, n_bins, qq) {
  pl1 <- ggplot(df, aes(x=rowMean, y=CV2, color=factor(Class))) +
    geom_point(aes(alpha = 0.5)) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("Mean expression, normalized counts") +
    ylab("CV^2") +
    stat_smooth(aes(group=factor(Class)), col=NA, method = "loess", size=1, se=TRUE) +
    stat_smooth(aes(color=factor(Class)), method = "loess", size=1, se=F) +
    theme(legend.justification=c(1,1), legend.position=c(0.99,0.99),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    scale_color_manual(values=c("#FC6621", "#105CFB"), name ="group") +
    scale_alpha(guide = 'none') +
    geom_vline(xintercept = qq[2:n_bins], linetype = "longdash", color="black")
  return(pl1)
}


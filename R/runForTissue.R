#' Generate plots for a given tissue
#'
#' @param name tissue name
#' @return
#' @export
#'
generatePlots <- function(name) {
  df <- process_to_counts(name)
  plot_CV2_vs_mean(df, name, plotMode = "Merged", modeBoxPlot = "Sizes", n_bins = 5)
}


# tissues <- c(
#   "Kidney_Cortex", "Heart_Left_Ventricle",
#   "Liver", "Lung", "Pancreas", "Stomach", "Small_Intestine", "Spleen"
# )
# for (fullName in tissues) {
#   generatePlots(fullName)
# }


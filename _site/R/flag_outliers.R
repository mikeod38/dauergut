#' flag_outliers
#'
#' Identify highly influential observations based on both standardized
#' residuals and cooks distance. Function is a wrapper for plot_outlier.resids
#' and plot_cooksD. 
#' @param lin.mod model fit, lm, glm, lmm, glmm
#' @param noplot optional to plot outlier observations, leave blank for no plot
#' @param threshold threshold for outlier flagging. Usually 4 (4*mean cooks D)
#' @param df data frame in which to flag outlier observations. Must be same as used in lin.mod
#' @return a list of observation indices based on the original data frame
#' @export
#' @examples df <- flag_outliers(lin.mod, df, threshold = 4, noplot = TRUE)
flag_outliers <- function(lin.mod, df, threshold, noplot) {
  outliers.q <- dauergut::plot_outlier.resids(lin.mod, noplot)
  outliers.cooksd <- dauergut::plot_cooksD(lin.mod, threshold, noplot)
  outliers <- intersect(outliers.q, outliers.cooksd)
  if(is.na(outliers[1])) {
    df$outlier.status <- FALSE
    newdata = df
  } else {
    df.outliers <- df[outliers,]
    df.outliers$outlier.status <- TRUE
    newdata <- left_join(df, df.outliers)
    newdata$outlier.status <- !is.na(newdata$outlier.status) #recode NA to FALSE
  }
  return(newdata)
}
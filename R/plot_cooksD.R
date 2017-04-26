#' plot_cooksD
#'
#' Identify highly influential observations based on cooks distance. 
#' Common outlier threshold is 4*mean of cooksD
#' @param lin.mod model fit, lm, glm, lmm, glmm
#' @param noplot optional to plot outlier observations, leave blank for no plot
#' @param threshold threshold for outlier flagging. Usually 4 (4*mean cooks D)
#' @return a list of observation indices based on the original data frame
#' @export
#' @examples outliers.cooksd <- plot_cooksD(fit, noplot = TRUE)
plot_cooksD<-function(lin.mod, threshold, noplot) {
  cooksd <- cooks.distance(lin.mod)
  if (missing(noplot)) {
    par(mfrow = c(1,1))
    plot(cooksd~lin.mod$fitted.values, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
    abline(h = threshold*mean(cooksd, na.rm=T), col="red")
    text(x=lin.mod$fitted.values-0.01*lin.mod$fitted.values, y=cooksd, 
         labels=ifelse(cooksd>threshold*mean(cooksd,na.rm=T),names(cooksd),""), col="red")  # add labels
    return(as.numeric(names(cooksd)[(cooksd > threshold*mean(cooksd, na.rm = T))]))
  } else {
    return(as.numeric(names(cooksd)[(cooksd >threshold*mean(cooksd, na.rm = T))]))
  }
}
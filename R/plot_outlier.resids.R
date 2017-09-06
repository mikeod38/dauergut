#' plot_outlier.resids
#'
#' identify outliers based on standardized residuals using MASS package. outliers are flagged if absolute standardized
#' residual is greater than 3 and is > 1.5*IQR for standardized residuals. 
#' @param lin.mod model fit, lm, glm, lmm, glmm
#' @param noplot optional to plot outlier observations, leave blank for no plot
#' @return a list of observation indices based on the original data frame
#' @export
#' @examples outliers.q <- plot_outlier.resids(fit, noplot = TRUE)
plot_outlier.resids <- function (lin.mod, noplot) {
  resid <- MASS::stdres(lin.mod)
  absresid<-abs(resid)
  if (missing(noplot)) {
    par(mfrow = c(1,2))
    plot(absresid~lin.mod$fitted.values, pch="*", cex=2, main="Outlier observations by standardized residuals")
    abline(h = 3, col="red")
    text(x=lin.mod$fitted.values-0.01*lin.mod$fitted.values,
         y=absresid, labels=ifelse(absresid > 3, names(absresid), ""),
         col ="red")
    plot(resid~lin.mod$fitted.values, main = "Outliers greater than 1.5*IQR for standardized residuals")
    abline(h = 1.5*IQR(resid), col = "blue", lty=2)
    abline(h = -1.5*IQR(resid), col = "blue", lty=2)
    text(x=lin.mod$fitted.values-0.01*lin.mod$fitted.values,
         y=resid,labels=ifelse(abs(resid) > 1.5*IQR(resid), names(resid), ""))
    return(as.numeric(names(absresid)[absresid > 1.5*IQR(resid)])) 
  } else {
    return(as.numeric(names(absresid)[absresid > 1.5*IQR(resid)]))
  }
}
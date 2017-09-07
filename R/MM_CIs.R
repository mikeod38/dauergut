#' MM_CIs
#'
#' generate Tukey-adjusted 95 pct confidence intervals for plotting on the
#' untransformed response scale for linear models
#' using lsmeans package
#' @param y input if model fit ie lm, glm, glmer etc...
#' @param predictor optional predictor variable for comparisons, defaults to "genotype"
#' @importFrom magrittr "%>%"
#' @export
#' @examples mixed <- MM_CIs(fit, predictor = "genotype")
MM_CIs<-function(y,predictor) {
  if(missing(predictor)) {
    predictor <- "genotype"
  } else {
  }
  strains <- y$xlevels
  y.rg<-lsmeans::ref.grid(y)
  y.lsm<-lsmeans::lsmeans(y, predictor)
  y.lsm.sum<-summary(y.lsm,level = .95, infer=TRUE, adjust = "Tukey")
  # mixed<-with(y.lsm.sum,data.frame(lsmean,lower.CL, upper.CL))
  colnames(y.lsm.sum)<-c(predictor,"mean", "SE", "df", "lower.CL", "upper.CL", "z.ratio", "p.value")
  mixed<-y.lsm.sum %>% dplyr::select(1, mean, lower.CL, upper.CL,SE)
  mixed$x.pos<-as.numeric(as.factor(mixed[[predictor]])) + 0.3
  return(mixed)
} 
#' MM_CI_trans_dunnett
#'
#' generate Dunnett-adjusted 95 pct confidence intervals for plotting on the
#' untransformed response scale for binomial or other transformed model data
#' using lsmeans package
#' @param y input if model fit ie lm, glm, glmer etc...
#' @param factor optional predictor variable for comparisons, defaults to "genotype"
#' @importFrom magrittr "%>%"
#' @export
#' @examples mixed <- MM_CI_trans_dunnett(fit, predictor = "genotype")
MM_CI_trans_dunnett<-function(y,predictor) {
  if(missing(predictor)) {
    predictor <- "genotype"
  } else {
  }
  strains <- levels(y@frame[[predictor]]) #grab levels of predictor factor from model
  y.rg<-lsmeans::ref.grid(y, type="response")
  y.lsm<-lsmeans::lsmeans(y.rg, predictor)
  y.lsm.sum<-summary(y.lsm,level = .95, infer=TRUE, adjust = "mvt", type="response")
  colnames(y.lsm.sum)<-c(predictor,"mean", "SE", "df", "lower.CL", "upper.CL", "z.ratio", "p.value")
  mixed<-y.lsm.sum %>% dplyr::select(1, mean, lower.CL, upper.CL,SE)
  mixed$x.pos<-as.numeric(as.factor(mixed[[predictor]])) + 0.3
  return(mixed)
}
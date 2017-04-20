#' tukey_contrasts
#'
#' generate post-hoc pairwise Tukey comparisons based on a model fit
#' results are presented on the original response scale
#' @param x model fit of type lm, lmer, glmer, glm etc...
#' @param factor primary predictor factor for comparisons ie genotype, temperature etc...
#' @return an lsmeans contrast list
#' @export
#' @examples contrasts <- tukey_contrasts(glmm.fit, "genotype")
tukey_contrasts<-function(x, factor) {
  #get post-hoc pairwise Tukey comparisons
  library(lsmeans)
  x.rg<-ref.grid(x, type = "response")
  contrasts<-x.rg %>% lsmeans(factor) %>% pairs(adjust = "mvt") %>% summary() %>% prange()
  return(contrasts)
}

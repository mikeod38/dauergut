#' tukey_contrasts
#'
#' generate post-hoc pairwise Tukey comparisons based on a model fit
#' results are presented on the original response scale
#' @param x model fit of type lm, lmer, glmer, glm etc...
#' @param factor primary predictor factor for comparisons ie genotype, temperature etc...
#' @return an lsmeans contrast list
#' @importFrom magrittr "%>%" 
#' @export
#' @examples contrasts <- tukey_contrasts(glmm.fit, "genotype")
tukey_contrasts<-function(x, factor) {
  #get post-hoc pairwise Tukey comparisons
  x.rg<-lsmeans::ref.grid(x, type = "response")
  if(length(strains) < 10) {
  contrasts<-x.rg %>% lsmeans::lsmeans(factor) %>% pairs(adjust = "mvt") %>% summary() %>% prange()
  } else {
    contrasts<-x.rg %>% lsmeans::lsmeans(factor) %>% pairs(adjust = "fdr") %>% summary() %>% prange()
  }
  return(contrasts)
}

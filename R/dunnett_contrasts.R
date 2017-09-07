#' dunnett_contrasts
#'
#' generate Dunnett-adjusted multiple comparisons
#' using lsmeans package
#' @param x input if model fit ie lm, glm, glmer etc...
#' @param ref.index index of reference level for multiple comparisons. You can find this
#' using names(<fit>$coefficients) for lm or glm obj or names(fixef(<fit>)) for mixed models
#' @param interaction optional term for interaction predictor if relevant. 
#' @param factor predictor variable for comparisons, usually a group ie genotype
#' @importFrom magrittr "%>%"
#' @export
#' @examples contrasts<-dunnett_contrasts(fit, ref.index = 1, factor = "genotype", interaction = "food")
dunnett_contrasts<-function(x, ref.index, factor, interaction) {
  #generates Dunnett contrasts for a given level of factor (use quotes for factor = ""), adjustment 
  x.rg <- lsmeans::ref.grid(x, type="response")
  contrasts.1<-x.rg %>% lsmeans::lsmeans(factor) %>% lsmeans::contrast("trt.vs.ctrl", ref = ref.index) %>% 
    summary(adjust = "mvt") %>%
    prange()
  if(missing(interaction)) {
    print("interaction term not indicated")
    return(contrasts.1)
  } else {
    contrasts.2<- x.rg %>% lsmeans::lsmeans(interaction, by = factor) %>% lsmeans::pairs(by = factor) %>%
      summary(adjust = "mvt", by=NULL) %>% prange()
    return(list(factor = contrasts.1, interaction = contrasts.2))
  }
}
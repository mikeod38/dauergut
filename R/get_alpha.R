#` get_alpha
#'
#'  summarize results of model simulations to get proportion p <0.05
#'  doesn't work if too few simulations and no p < 0.05 for one of the models
#' 
#' 
#' @param df simulation data frame from sim_dauer, sim_dauer_unbal, or sim_dauer_biased. Pre-processed
#' using add.setts.attr to add attributes
#' @param nsim
#' 
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' 
#' @export
#' @examples balanced.mod %<>% get_alpha()

get_alpha <- function(df) {
  # for t
  t1 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype2 < 0.05)
  t2 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype3 < 0.05)
  
  # for anova
  lm1 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype2 < 0.05 & Fp < 0.05)
  lm2 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype3 < 0.05 & Fp < 0.05)
  
  # for glmm
  glmm1 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype2 < 0.05 & Chisq.p < 0.05)
  glmm2 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype3 < 0.05 & Chisq.p < 0.05)
  
  #for stan
  stan1 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype2 < 0.05)
  stan2 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype3 < 0.05)
  
  nsim <- max(df$dataset)
  
  alpha = list(t1 = t1,t2 = t2,lm1 = lm1, lm2 = lm2, glmm1 = glmm1, glmm2 = glmm2, stan1 = stan1, stan2 = stan2)
  
  prop_sig <- function(df) {
    # number sig in contingency table = [2,2]
    df[2,2]/nsim
  }
  alpha = c(lapply(alpha, prop_sig), nsim = nsim)
  return(alpha)
}

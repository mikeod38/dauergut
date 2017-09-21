#' add.setts.attr
#'
#'add attributes to simualtion based on input settigns
#' 
#' 
#' @param simulation simulation from sim_dauer, sim_dauer_unbal, or sim_dauer_biased
#' requires a 'settings' list in the environment
#' @param settings list of settings from simluations
#' 
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr "%>%"
#' 
#' @export
#' @examples simulation %<>% add.setts.attr
add.setts.attr <- function(simulation, settings) {
  simulation %<>% mutate(dataset = rep(1:(nrow(.)/length(levels(simulation$model))),each = length(levels(simulation$model))),
                         genotype2=as.numeric(as.character(genotype2)),
                         genotype3=as.numeric(as.character(genotype3)),
                         Fp=as.numeric(as.character(Fp)),
                         `Chisq.p` = as.numeric(as.character(`Chisq.p`))
  ) %>% `attr<-`('settings', unlist(settings))
}

#' par.replicate
#'
#'  wrapper fo parSapply to mimic replicate in base R for parallel
#' 
#' 
#' @param cl cluster connections from parallel::makeCluster
#' @param n number of replicates, or number of simulations
#' @param expr expression or function to use in replicate
#' @param simplify wrapper to simplify within parSapply, usually FALSE
#' 
#' @export
#' @examples library(parallel) 
#' cl <- makeCluster(2) 
#' clusterEvalQ(cl,  { library(MASS)
#' library(magrittr)
#' library(dplyr)
#' library(lme4)
#' library(rstan)
#' library(rstanarm)
#' })
#' clusterExport(cl=cl, varlist=c("sim_dauer_3G_stan", "settings"))
#' simulation <- do.call( rbind, par.replicate(cl,n=1000, 
#' sim_dauer_3G_stan(c(settings,
#' do.plot = FALSE, 
#' do.stan = TRUE)), 
#' simplify=FALSE ))
#' 
par.replicate <- function(cl, n, expr, simplify = FALSE) {
  parSapply(cl, integer(n), eval.parent(substitute(function(...) expr)), 
            simplify = simplify)
}
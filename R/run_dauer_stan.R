#' run_dauer_stan
#'
#' Function runs a stan glmm for dauer data
#' 
#' @param df input dataset. Requires a "genotype" column. See type for data column types
#' @param type "dauer" (default),"dauer-grouped". For dauer data, needs to have raw counts, with "dauer" column and "n" column.
#' For grouped data, must have group.id column - usually interaction(genotype,condition)
#' 
#' @export
#' @examples  df %>% run_dauer_stan(parameters)

run_dauer_stan <- function(df,type,parameters) {

  rstan::rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) # Run on multiple cores
  
  if(missing(parameters)) {
    parameters = list(chains = 3, cores =4, seed = 2000,iter=6000,
    control = list(adapt_delta=0.99))
  } else {
    parameters = parameters
  }
  
if(missing(type)) {
  mod <- rstanarm::stan_glmer(data=df,
  formula =  cbind(dauer, (n-dauer)) ~ genotype + (1|day) + (1|strainDate) + (1|plateID),
  family = binomial(link="logit"),
  chains = parameters$chains,
  cores = parameters$cores,
  seed = parameters$seed,
  iter = parameters$iter,
  control = list(adapt_delta=0.99))
} else {
  if(type == "dauer-grouped") {
    mod <- rstanarm::stan_glmer(formula = cbind(dauer, (n-dauer)) ~ 0 + group.id + (1|day) + (1|strainDate) + (1|plateID),
               data=df,
               family = binomial(link="logit"),
               chains = parameters$chains, 
               cores = parameters$cores, 
               seed = parameters$seed,
               iter=parameters$iter,
               control = list(adapt_delta=0.99))
  } else {
    print("invalid type, write full model")
  }
  return(mod)
}
}
  
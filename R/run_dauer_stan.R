#' run_dauer_stan
#'
#' Function runs a stan glmm for dauer data
#' 
#' @param df input dataset. Requires a "genotype" column. See type for data column types
#' @param type "dauer" (default),"dauer-grouped". For dauer data, needs to have raw counts, with "dauer" column and "n" column.
#' For grouped data, must have group.id column - usually interaction(genotype,condition)
#' @param group for grouped data, specify the name of the grouping factor
#' @param intercept optional argument for intercept in the model. use for direct estimation of the reference group. default is TRUE
#' 
#' @export
#' @examples  df %>% run_dauer_stan(parameters)
#' @examples  run_dauer_stan(df, type = "dauer-grouped", group = "food", intercept = "false")

run_dauer_stan <- function(df,type,parameters,group,intercept) {

  rstan::rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) # Run on multiple cores
  
  if(missing(parameters)) {
    parameters = list(chains = 3, cores =4, seed = 2000,iter=6000,
    control = list(adapt_delta=0.99))
  } else {
    parameters = parameters
  }
  
  
if(missing(group)) {
  formula =as.formula("cbind(dauer, (n-dauer)) ~ genotype + (1|day) + (1|plateID) + (1|strainDate)")
} else {
  if(type == "dauer-grouped" & missing(intercept))  {
    formula = as.formula(paste("cbind(dauer, (n-dauer)) ~",
                               "group.id + (1|day) + (1|plateID)", "+ (1|strainDate)"))
  } else {
    if(type == "dauer-grouped" & intercept == FALSE) {
      formula = as.formula(paste("cbind(dauer, (n-dauer)) ~",
                                 "0 + group.id + (1|day) + (1|plateID)", "+ (1|strainDate)"))
    } else {
      print("invalid type, write full model")
    }
  }
}
  
  

  mod <- rstanarm::stan_glmer(data=df,
  formula =  formula, 
  family = binomial(link="logit"),
  chains = parameters$chains,
  cores = parameters$cores,
  seed = parameters$seed,
  iter = parameters$iter,
  control = list(adapt_delta=0.99))

  return(mod)
  
  }


  
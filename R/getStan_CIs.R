#' getStan_CIs
#'
#' Function pulls 50 and 95% confidence intervals from a StanGLMM model using genotype as a parameter/predictor
#' @param x Stan glmm model
#' @param type type of data, options are NA (linear model), "log" (Log10), "roam" (Poisson), "dauer" (logit)
#' @param group optional argument for group-level effects, i.e. "food" or "temp". Must match column name for condition
#' @param base optional argument for log models - base of log function
#' @param compare_gt optional argument for group-level effect - if true, genotype is compared at each level of condition, if blank,
#' then condition effects are returned.
#' @param intercept optional parameter specifying whether an intercept was included in the model formula - see run_dauer_stan.R. 
#' If blank, assumes there is an intercept in model.
#' 
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples getStan_CIs(stan.glmm, type = "log", group = "food", base = 10, compare_gt = TRUE)
#' 
getStan_CIs <- function(x,type,group,base,compare_gt,intercept) {
  
  # get estimated parameters from the model:
  if(missing(group)) {
    strains = levels(x$data$genotype) # for non-grouped data
    data <- summary(x)[1:length(strains),c(1,4,8,5,7)] + summary(x)[1,1]
    data[1,] <- data[1,] - summary(x)[1,1] #same for  all non-grouped data
  } else {
    if(missing(intercept)) {
    group.id <- levels(x$data$group.id) # same for all grouped data
    my_group <- x$data[,group]
    data <- summary(x)[1:length(group.id),c(1,4,8,5,7)] + summary(x)[1,1]
    data[1,] <- data[1,] - summary(x)[1,1]
    } else {
      group.id <- levels(x$data$group.id) # no intercept so direct estimates of cred intervals. 
      my_group <- x$data[,group]
      data <- summary(x)[1:length(group.id),c(1,4,8,5,7)]
    }
  }

  # return either transformed or non-transformed data
  if(missing(type)) {
    data %<>% data.frame()
  } else { 
    if(type == "log") {
      data <- base^data %>% data.frame()
    } else {
      if(type == "roam") {
        data %<>% exp() %>% data.frame()
      } else {  # for type == "dauer"
        data %<>% brms::inv_logit_scaled() %>% data.frame() 
      }
    }
  }
  
  # format output:
  
  if(missing(group)) {
    data %<>% mutate(genotype = strains) %>% `colnames<-`(c("mean", "lower.CL", "upper.CL", "lower.25", "upper.75", "genotype")) %>%
      mutate(x.pos = seq(1:length(strains)) + 0.3)
  } else {
    data %<>% mutate(genotype = rep(strains, 2),
                     group = rep(levels(my_group), each = length(strains)),
                     group.id = group.id) %>% `colnames<-`(c("mean", "lower.CL", "upper.CL", "lower.25", "upper.75", "genotype", group, "group.id"))
    if(missing(compare_gt)) {
      data %<>% mutate(x.pos = rep(1:2, each = length(strains)) + 0.3)
    } else {
      data %<>% mutate(x.pos = rep(1:2, length(strains)) + 0.3)
    }
  }
  return(data)
}
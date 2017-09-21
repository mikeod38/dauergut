#' sim_dauer
#'
#' simulate balanced dauer assay data using 3 groups of strains
#' I, A and B. Simulation uses day-to-day and plate-to-plate variance. 
#' Can be used to simulate any binomial data with hierarchical clusters.
#' 
#' 
#' @param settings input list of settings for the simulation. Must be a list.
#' settings$I    #population control intercept (in logit)
#' settings$nP  # number of plates
#' settings$nD  # number of days
#' settings$sP  # plate to plate sd
#' settings$sD  # day to day sd
#' settings$sG  # genotype sd due to culture history (logit)
#' settings$k   # number animal per plate)
#' settings$A    # population A intercept (expt)
#' settings$B    # pop B intercept
#' 
#' @export
#' @examples settings <- list(settings <- list(
#' I = 0 #population control intercept (in logit). 0 = p(0.5)
#' ,nP = 6 # number of plates
#' ,nD = 3 # number of days
#' ,sP = 0.1 # plate to plate variance (0.3)
#' ,sD = 0.5 # day to day variance (0.2)
#' ,sG = 0.5 # genotype variance due to culture history (logit) (0.2)
#' ,k = 60 # number animal per plate)
#' ,A = 0 # population A intercept (expt - genotype2)
#' ,B = 0 #pop B intercept (expt - genotype2)
#' ))
#' 
#' 
#' sim_dauer(settings = c(settings, do.plot = TRUE))
#' 
#' simulation <- sim_dauer(settings = c(settings, do.plot = FALSE, do.stan = TRUE))

sim_dauer<-function(settings) {
  # get settings
  I = settings$I    #population control intercept (in logit)
  nP = settings$nP  # number of plates
  nD = settings$nD  # number of days
  sP <- settings$sP  # plate to plate sd
  sD = settings$sD  # day to day sd
  sG = settings$sG  # genotype sd due to culture history (logit)
  k  = settings$k   # number animal per plate)
  A = settings$A    # population A intercept (expt)
  B = settings$B    # pop B intercept
  do.plot = settings$do.plot # plot for vis inpection (no model fits)
  do.stan = settings$do.stan # fit stan_glmer models
  
  ############# model functions #######################
  lm.sim <-   function(df) {
    modsum <- df %>% lm(formula = p~genotype) %>% summary()
    genotype2 <- as.numeric(modsum$coefficients[,4][2])
    genotype3 <- as.numeric(modsum$coefficients[,4][3])
    Fp <- as.numeric(1-pf(modsum$fstatistic[1],modsum$fstatistic[2], modsum$fstatistic[3]))
    Chisq.p = NA
    model <- "anova"
    p.val <- data.frame(cbind(model, genotype2, genotype3, Fp, Chisq.p))
    return(p.val)
  }
  
  t.sim <- function(df) {
    genotype2 = data %>% dplyr::filter(genotype != "3") %$% t.test(p~genotype)$p.value
    genotype3 = data %>% dplyr::filter(genotype != "2") %$% t.test(p~genotype)$p.value
    model = "t"
    Fp = NA
    Chisq.p = NA
    p.val <- data.frame(cbind(model,genotype2, genotype3, Fp, Chisq.p))
    return(p.val)
  }
  
  glmm.sim <- function(df) {
    mod = data %>%
      lme4::glmer(formula = cbind(y, (k-y)) ~ genotype + (1|day/strainDate/plateID),
                  family = binomial, control=glmerControl(optimizer="bobyqa"))
    nullmod = data %>% lme4::glmer(formula = cbind(y, (k-y)) ~ 1 + (1|day/strainDate/plateID),
                                   family = binomial, control=glmerControl(optimizer="bobyqa"))
    modsum <- mod %>% summary()
    genotype2 <- as.numeric(modsum$coefficients[,4][2])
    genotype3 <- as.numeric(modsum$coefficients[,4][3])
    model <- "glmm"
    compmod <- anova(nullmod, mod)
    Fp = NA
    Chisq.p <- compmod$`Pr(>Chisq)`[2]
    p.val <- data.frame(cbind(model, genotype2, genotype3, Fp, Chisq.p))
    return(p.val)
  }
  
  stan.sim <- function(df) {
    library (rstan)
    rstan_options (auto_write=TRUE)
    options (mc.cores=parallel::detectCores ()) # Run on multiple cores
    # run stan mod with default priors
    mod <- stan_glmer( cbind(y, k-y) ~ genotype + (1|day) + (1|strainDate) + (1|plateID),
                       data=data,
                       family = binomial(link="logit"),
                       chains = 3, cores =4, seed = 2000,
                       control = list(adapt_delta=0.99)
    )
    model = "stan"
    # get posterior 95% cred interval, test if it contains 0 (abs(sum) != sum(abs))
    mod.pp <- posterior_interval(mod, prob = 0.95, pars = c("genotype2", "genotype3"))
    # will give TRUE if 95% CI contains 0
    genotype2 <- as.numeric(abs(mod.pp[1,1]) + abs(mod.pp[1,2]) != abs(mod.pp[1,1] + mod.pp[1,2]))
    genotype3 <- as.numeric(abs(mod.pp[2,1]) + abs(mod.pp[2,2]) != abs(mod.pp[2,1] + mod.pp[2,2]))
    Fp = NA
    Chisq.p = NA
    p.val <- data.frame(cbind(model, genotype2, genotype3, Fp, Chisq.p))
    return(p.val)
    #return(mod)
  }
  # 
  #   
  #   ############### generate simulated data #############
  gen.dauer.data <- function(...) {
    # random effects with mean 0 and var = sP,sD or sG N
    RE.p.I = as.numeric(rnorm(nP, 0, sd = sP))
    RE.p.A = as.numeric(rnorm(nP, 0, sd = sP))
    RE.p.B = as.numeric(rnorm(nP, 0, sd = sP))
    RE.GP.I = as.numeric(rnorm(nD, 0, sd = sG))
    RE.GP.A = as.numeric(rnorm(nD, 0, sd = sG))
    RE.GP.B = as.numeric(rnorm(nD, 0, sd = sG))
    RE.d = as.numeric(rnorm(nD, 0, sd = sD))
    
    
    day = (seq(1:nD))
    plate = seq(1:nP)
    
    # data for three groups - balanced
    data.I <- cbind(genotype = 1,
                    plate = plate,
                    mean = I,
                    k = rpois(nP, 60),
                    RE.p = RE.p.I,
                    day = rep(day, each = nP/nD),
                    RE.d = rep(RE.d, each = nP/nD),
                    RE.GP = rep(RE.GP.I, each = nP/nD),
                    y = NA) %>% data.frame() %>%
      dplyr::mutate(y=rbinom(nP,k,boot::inv.logit(RE.p + RE.d + RE.GP + mean)))
    
    data.A <- cbind(genotype = 2,
                    plate = plate,
                    mean = A,
                    k = rpois(nP, 60),
                    RE.p = RE.p.A,
                    day = rep(day, each = nP/nD),
                    RE.d = rep(RE.d, each = nP/nD),
                    RE.GP = rep(RE.GP.A, each = nP/nD),
                    y = NA) %>% data.frame() %>%
      dplyr::mutate(y=rbinom(nP,k,boot::inv.logit(RE.p + RE.d + RE.GP + mean)))
    
    data.B <- cbind(genotype = 3,
                    plate = plate,
                    k = rpois(nP, 60),
                    mean = B,
                    RE.p = RE.p.B,
                    day = rep(day, each = nP/nD),
                    RE.d = rep(RE.d, each = nP/nD),
                    RE.GP = rep(RE.GP.B, each = nP/nD),
                    y = NA) %>% data.frame() %>%
      dplyr::mutate(y=rbinom(nP,k,boot::inv.logit(RE.p + RE.d + RE.GP + mean)))
    
    data <- rbind(data.I, data.A, data.B) %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    strainDate = interaction(genotype,day),
                    plateID = interaction(genotype,plate),
                    p = y/k)
    return(data)
  }
  
  data <- gen.dauer.data()
  
  # optional plot (use only for single simulation inspection)
  if(do.plot) {
    p<-data %>% ggplot(aes(x=genotype, y=p)) +
      geom_boxplot() +
      geom_point(aes(x=genotype, colour = factor(day)))
    return(p)
  } else {
    if(do.stan) {
      lm <- lm.sim(data)
      t <- t.sim(data)
      glmm <- glmm.sim(data)
      stan <- stan.sim(data)
      
      p.val <- rbind(lm, t, glmm, stan)
      return(p.val)
    } else {
      lm <- lm.sim(data)
      t <- t.sim(data)
      glmm <- glmm.sim(data)
      
      p.val <- rbind(lm, t, glmm)
      return(p.val)
    }
  }
}
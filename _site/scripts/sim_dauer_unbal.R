########## simulation function ####################
sim_dauer_unbal<-function(settings) {
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
  do.stan = settings$do.stan # fit stan_glmer 

  
  ############### generate simulated data #############
  
  day = (seq(1:nD))
  plate = seq(1:(nP))
  
  missing.days.1 <- sample(day,2)
  missing.days.2 <- sample(day[!day %in% missing.days.1],2)
  
  gen.dauer.data <- function(...) {
    # random effects with mean 0 and var = sP,sD or sG N
    RE.p.I = as.numeric(rnorm(nP, 0, sd = sP))
    RE.p.A = as.numeric(rnorm(nP, 0, sd = sP))
    RE.p.B = as.numeric(rnorm(nP, 0, sd = sP))
    RE.GP.I = as.numeric(rnorm(nD, 0, sd = sG))
    RE.GP.A = as.numeric(rnorm(nD, 0, sd = sG))
    RE.GP.B = as.numeric(rnorm(nD, 0, sd = sG))
    RE.d = as.numeric(rnorm(nD, 0, sd = sD))
    
    # data for three groups - balanced
    data.I <- cbind(genotype = 1,
                    plate = plate,
                    mean = I,
                    k = rpois(nP, k),
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
                    y = NA) %>% data.frame() %>% dplyr::filter(!day %in% missing.days.1) %>%
      dplyr::mutate(y=rbinom(nP-(length(missing.days.1)*2),k,boot::inv.logit(RE.p + RE.d + RE.GP + mean)))
    
    data.B <- cbind(genotype = 3,
                    plate = plate,
                    k = rpois(nP, 60),
                    mean = B,
                    RE.p = RE.p.B,
                    day = rep(day, each = nP/nD),
                    RE.d = rep(RE.d, each = nP/nD),
                    RE.GP = rep(RE.GP.B, each = nP/nD),
                    y = NA) %>% data.frame() %>% dplyr::filter(!day %in% missing.days.2) %>%
      dplyr::mutate(y=rbinom(nP-(length(missing.days.2)*2),k,boot::inv.logit(RE.p + RE.d + RE.GP + mean)))
    
    data <- rbind(data.I, data.A, data.B) %>%
      dplyr::mutate(genotype = as.factor(genotype),
                    strainDate = interaction(genotype,day),
                    plateID = interaction(genotype,plate),
                    p = y/k)
    return(data)
  }
  
  data <- gen.dauer.data()
  
  
    
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
    genotype2 = data %>% dplyr::filter(genotype != "3" & !day %in% missing.days.1) %$% t.test(p~genotype)$p.value
    genotype3 = data %>% dplyr::filter(genotype != "2" & !day %in% missing.days.2) %$% t.test(p~genotype)$p.value
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



####### run test simulations ############################
settings <- list(I = -2.2 #population control intercept (in logit). 0 = p(0.5)
                 ,nP = 10 # number of plates in total
                 ,nD = 5 # number of days (some will be missing)
                 ,sP = 0.3 # plate to plate variance
                 ,sD = 0.2 # day to day variance
                 ,sG = 0.2 # genotype variance due to culture history (logit)
                 ,k = 60 # number animal per plate)
                 ,A = -2.2  # population A intercept (expt - genotype2)
                 ,B = -2.2 #pop B intercept (expt - genotype3)
)

#check data distributions with sample simulations
(p <- sim_dauer_unbal(c(settings, do.plot = TRUE, do.stan = FALSE)))
(simulation <- sim_dauer_unbal(c(settings, do.plot = FALSE, do.stan = FALSE)))

############ for parallel sampling below below ###########
library(parallel)
cl <- makeCluster(6) ### stan uses 3 cores(3 chains) so total = cl * 3
# 2 clusters (6 cores) = ~4sec/sim
clusterEvalQ(cl,  { library(MASS)
  library(magrittr)
  library(dplyr)
  library(lme4)
  library(rstan)
  library(rstanarm)
})
settings$do.plot = FALSE
clusterExport(cl=cl, varlist=c("sim_dauer_unbal", "settings"))

# wrapper fo parSapply to mimic replicate in base R
par.replicate <- function(cl, n, expr, simplify = FALSE) {
  parSapply(cl, integer(n), eval.parent(substitute(function(...) expr)), 
            simplify = simplify)
}

# sampling - n simulations ~ 1hr using 6 cores
system.time(simulation <- do.call( rbind, par.replicate(cl,n=1000, sim_dauer_unbal(c(settings, do.plot = FALSE, do.stan = FALSE)), simplify=FALSE )))
# output is a list of p.values (and/or binary cutoff with alpha < 0.05)
# append attributes format to numeric
simulation %<>% mutate(dataset = rep(1:(nrow(.)/length(levels(sim$model))),each = length(levels(sim$model))),
                genotype2=as.numeric(as.character(genotype2)),
                genotype3=as.numeric(as.character(genotype3)),
                Fp=as.numeric(as.character(Fp)),
                `Chisq.p` = as.numeric(as.character(`Chisq.p`))
) %>% `attr<-`('settings', unlist(settings[1:9]))



#analysis of alpha frequency                                                             
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
  
  alpha = list(t1 = t1,t2 = t2,lm1 = lm1, lm2 = lm2, glmm1 = glmm1, glmm2 = glmm2, stan1 = stan1, stan2 = stan2)
  
  prop_sig <- function(df) {
    # number sig in contingency table = [2,2]
    df[2,2]/1000
  }
  alpha = lapply(alpha, prop_sig)
  return(alpha)
}

#### compile a list of model p-value outputs from simulations #####
#make list of simulation outcomes (must keep in workspace)
unbalanced.mods = mget(ls(pattern = "simulation"))
out <- lapply(unbalanced.mods,function(x) {x %>% get_alpha %>% unlist})
# save output to ASCI
dput(out, "data/sim_data/model_sim_alphas")
dput(balanced.mods, "data/sim_data/model_simdata")
#make table of results
#for H0: (G1=G2=G3) = TRUE
mget(ls(pattern = "b0")) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
  t() %>% `row.names<-`(c("mean = 0.5", "mean = 0.5, highVar", "mean = 0.1", "mean = 0.1, highVar")) %>%
  `colnames<-`(c("t1", "t2", "lm1", "lm2", "glmm1", "glmm2", "stan.glmm1", "stan.glmm2")) %>% kable()

#for H0: (G1=G2=G3) = FALSE, H1: G3 > G1 = G2
mget(ls(pattern = "b1")) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
  t() %>% `row.names<-`(c("mean = 0.5", "mean = 0.5, highVar", "mean = 0.1", "mean = 0.1, highVar")) %>%
  `colnames<-`(c("t1", "t2", "lm1", "lm2", "glmm1", "glmm2", "stan.glmm1", "stan.glmm2")) %>% kable()

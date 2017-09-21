
####### run test simulations ############################
settings <- list(I = 0 #population control intercept (in logit). 0 = p(0.5)
                 ,nP = 6 # number of plates
                 ,nD = 3 # number of days
                 ,sP = 0.1 # plate to plate variance (0.3)
                 ,sD = 0.5 # day to day variance (0.2)
                 ,sG = 0.5 # genotype variance due to culture history (logit) (0.2)
                 ,k = 60 # number animal per plate)
                 ,A = 0 # population A intercept (expt - genotype2)
                 ,B = 0 #pop B intercept (expt - genotype2)
)

#check data distributions with sample simulations
(p <- sim_dauer(c(settings, do.plot = TRUE)))
(sim <-sim_dauer(c(settings, do.plot = FALSE, do.stan = TRUE)))

############ for parallel sampling below below ###########
library(parallel)
cl <- makeCluster(2) ### stan uses 3 cores(3 chains) so total = cl * 3
# 2 clusters (6 cores) = ~4sec/sim
clusterEvalQ(cl,  { library(MASS)
library(magrittr)
library(dplyr)
library(lme4)
library(rstan)
library(rstanarm)
  library(dauergut)
  })
settings$do.plot = FALSE
clusterExport(cl=cl, varlist=c("sim_dauer", "settings"))


##### sampling ~ 700 simulations ~ 1hr using 6 cores ######
simulation <- do.call( rbind, 
                       par.replicate(cl,n=100, # number of sims here 
                                     sim_dauer(c(settings, # simulation fxn goes here
                                                 do.plot = FALSE,
                                                 do.stan = TRUE)),
                                     simplify=FALSE )) %>% add.setts.attr(settings = c(settings, model = "balanced"))


# output is a list of p.values (and/or binary cutoff with alpha < 0.05)

filename <- paste(names(attributes(simulation)$settings[c(1,4:6,8,9,11)]),
      attributes(simulation)$settings[c(1,4:6,8,9,11)], sep = "", collapse = "_")

dput(simulation, file.path(pathname, "extdata/sim_data/", paste(filename,"data", sep = ".",collapse = ".")))

# visually inspect results
(out <- simulation %>% get_alpha() %>% unlist)
dput(out, file.path(pathname, "extdata/sim_data/", paste(filename,"alphas", sep = ".",collapse = ".")))

# #### compile a list of model p-value outputs from simulations #####
# #make list of simulation outcomes
# simulated.mods = mget(ls(pattern = "1000x"))
# out <- lapply(simulated.mods,function(x) {x %>% get_alpha %>% unlist})
# out
# # save output to ASCI
# dput(out, "data/sim_data/model_sim_alphas")
# dput(simulated.mods, "data/sim_data/model_simdata")
# 
# 
# #make table of results
# #for H0: (G1=G2=G3) = TRUE
# mget(ls(pattern = "b0")) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
#   t() %>% 
#   `row.names<-`(c("mean = 0.5", "mean = 0.5, biased", "mean = 0.5, unbalanced", "mean = 0.5, highVar",
#                   "mean = 0.1", "mean = 0.1, biased", "mean = 0.1, unbalanced", "mean = 0.1, highVar")) %>%
#   `colnames<-`(c("t.1", "t.2", "lm.1", "lm.2", "glmm.1", "glmm.2", "stan.glmm.1", "stan.glmm.2")) %>% kable(caption = "for H0: (G1=G2=G3) = TRUE")
# 
# #for H0: (G1=G2=G3) = FALSE, H1: G3 > G1 = G2
# mget(ls(pattern = "b1")) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
#   t() %>% 
#   `row.names<-`(c("mean = 0.5", "mean = 0.5, biased", "mean = 0.5, unbalanced", "mean = 0.5, highVar",
#                   "mean = 0.1", "mean = 0.1, biased", "mean = 0.1, unbalanced", "mean = 0.1, highVar")) %>%
#   `colnames<-`(c("t.1", "t.2", "lm.1", "lm.2", "glmm.1", "glmm.2", "stan.glmm.1", "stan.glmm.2")) %>% kable(caption = "H0: (G1=G2=G3) = FALSE, H1: G3 > G1 = G2")




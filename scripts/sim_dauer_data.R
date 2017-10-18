
####### run test simulations ############################
settings <- list(I = 0 #population control intercept (in logit). 0 = p(0.5)
                 ,nP = 12 # number of plates
                 ,nD = 6 # number of days
                 ,sP = 0.1 # plate to plate variance (0.3)
                 ,sD = 0.5 # day to day variance (0.2)
                 ,sG = 0.5 # genotype variance due to culture history (logit) (0.2)
                 ,k = 60 # number animal per plate)
                 ,A = 0 # population A intercept (expt - genotype2)
                 ,B = 0 #pop B intercept (expt - genotype2)
)

#check data distributions with sample simulations
(p <- sim_dauer_biased(c(settings, do.plot = TRUE)))
(sim <-sim_dauer_biased(c(settings, do.plot = FALSE, do.stan = TRUE)))

############ for parallel sampling below below ###########
library(parallel)
cl <- makeCluster(2) ### stan uses 3 cores(3 chains) so total = cl * 3
# 2 clusters (6 cores) = ~4sec/sim
clusterEvalQ(cl,  { library(MASS)
library(magrittr)
library(dplyr)
library(lme4)
library(lmerTest)
library(lsmeans)
library(rstan)
library(rstanarm)
library(dauergut)
  })
clusterExport(cl=cl, varlist=c("sim_dauer", "sim_dauer_unbal", 
                               "sim_dauer_biased", "settings", "par.replicate",
                               "add.setts.attr"))


##### sampling ~ 700 simulations ~ 1hr using 6 cores ######
simulation <- do.call( rbind, 
                       par.replicate(cl,n=1000, # number of sims here 
                                     sim_dauer_biased(c(settings, # simulation fxn goes here
                                                 do.plot = FALSE,
                                                 do.stan = FALSE)),
                                     simplify=FALSE )) %>% 
  add.setts.attr(settings = c(settings, model = "biased")) # change to "unbal" or "biased" dep on fxn

stopCluster(cl)

#### output is a list of p.values (and/or binary cutoff with alpha < 0.05) ####
filename <- paste(names(attributes(simulation)$settings[c(1,4:6,8,9,10)]),
      attributes(simulation)$settings[c(1,4:6,8,9,10)], sep = "", collapse = "_")
dput(simulation, file.path(pathname, "extdata/sim_data/", paste(filename,"data", sep = ".",collapse = ".")))
#for multiple simulations, rename object with filename
assign(filename,simulation)

# visually inspect results
(out <- simulation %>% get_alpha() %>% unlist)
dput(out, file.path(pathname, "extdata/sim_data/", paste(filename,"alphas", sep = ".",collapse = ".")))


#### compile a list of model p-value outputs from simulations #####
#make list of simulation outcomes (must keep in workspace)
balanced.mods = mget(ls(pattern = "balanced"))
out <- lapply(balanced.mods,function(x) {x %>% get_alpha %>% unlist})
# save output to ASCI
dput(out, "data/sim_data/model_sim_alphas")
dput(balanced.mods, "data/sim_data/model_simdata")
#make table of subset of results
#for H0: (G1=G2=G3=0) = TRUE
mget(ls(pattern = "B0") %>% grep(pattern = "I0", value = TRUE)) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
  t() %>% `row.names<-`(c("mean = 0.5", "mean = 0.5, highVar", "mean = 0.1", "mean = 0.1, highVar")) %>%
  `colnames<-`(c("t1", "t2", "lm1", "lm2", "glmm1", "glmm2", "stan.glmm1", "stan.glmm2")) %>% kable()

#for H0: (G1=G2=G2=0) = FALSE, H1: G3 > G1 = G2, (G3 = G1+1)
mget(ls(pattern = "B1") %>% grep(pattern = "I0", value = TRUE)) %>% sapply(.,function(x) {x %>% get_alpha %>% unlist}) %>% 
  t() %>% `row.names<-`(c("mean = 0.5", "mean = 0.5, highVar", "mean = 0.1", "mean = 0.1, highVar")) %>%
  `colnames<-`(c("t1", "t2", "lm1", "lm2", "glmm1", "glmm2", "stan.glmm1", "stan.glmm2")) %>% kable()


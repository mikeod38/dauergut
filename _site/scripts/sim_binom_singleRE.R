library(MASS)
library(lme4)

generate_data = function(
  n # number of units/plates
  , k # number of trials within each condition within each unit
  , I # population intercept
  , vU # across-units variance of intercepts
  , A # population A effect
  #, vA # across-units variance of A effects
  #, rIA # across-units correlation between intercepts and A effects
){
  # uncorrelated random effects
  Sigma = vU
  Sigma = matrix(Sigma,1,1)
  means = mvrnorm(n,c(I),Sigma) # should code as normal mean 0 
  temp = expand.grid(A=c('a1','a2','a3'),value=0)
  temp$A = factor(temp$A)
  contrasts(temp$A) = contr.sum
  from_terms = terms(value~A)
  mm = model.matrix(from_terms,temp)
  data = expand.grid(A=c('a1','a2', 'a3'),unit=1:n,trial=1:k)
  for(i in 1:n){
    data$value[data$unit==i] = rbinom(k*3,1,plogis(as.numeric(mm %*% means[i,])))
  }
  data$unit = factor(data$unit)
  data$A = factor(data$A)
  contrasts(data$A) = contr.sum
  return(data)
}


this_data = generate_data(
  n = 6 # number of units
  , k = 60 # number of trials within each condition within each unit
  , I = 0.5 # population intercept
  , vU = .1 # across-units variance of intercepts
  , A = 0 # population A effect
  #, vA = .2 # across-units variance of A effects
  #, rIA = .6 # across-units correlation between intercepts and A effects
)

fit = lmer(
  data = this_data
  , formula = value ~ (1+A|unit) + A
  , family = binomial
)
print(fit)


this_data %>% group_by(A) %>% summarise(sum(value)/n())

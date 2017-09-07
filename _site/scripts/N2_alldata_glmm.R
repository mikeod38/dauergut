strains = "N2"
N2_alldata <- read.csv(file.path(pathname, "data", "N2_alldata.csv")) %>% format_dauer()

N2_alldata %>% ggplot(aes(x = month,y = pct, group = month)) + geom_boxplot() +
  geom_jitter(aes(colour = month), width = 0.05) + guides(legend=NULL)


N2.glmm.day <- glmer(data = N2_alldata, cbind(dauer,non) ~ 1 + (1|day), family=binomial)
N2.glmm.plate <- glmer(data = N2_alldata, cbind(dauer,non) ~ 1 + (1|plateID), family=binomial)
N2.glmm.day.plate <- update(N2.glmm.plate, formula = .~. + (1|day))
N2.month.n <- update(N2.glmm.day.plate, formula = .~. + n + factor(year), glmerControl(optimizer = "bobyqa"))
summary(N2.glmm.day.plate)

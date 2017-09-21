### script to analyze laarge dataset of N2 dauer formation
### execute from package dir with source(file.path(pathname, "scripts/N2_alldata_glmm.R"))
pathname = getwd()
strains = "N2"
N2_alldata <- read.csv(file.path(pathname, "extdata", "N2_alldata.csv")) %>% format_dauer(p.dauer = "non")

print(
  N2_alldata %>% ggplot(aes(x = genotype,y = pct)) + geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.025) + guides(legend=NULL)
)


N2.glmm.day <- glmer(data = N2_alldata, cbind(dauer,non) ~ 1 + (1|day), family=binomial)
N2.glmm.plate <- glmer(data = N2_alldata, cbind(dauer,non) ~ 1 + (1|plateID), family=binomial)
N2.glmm.day.plate <- update(N2.glmm.plate, formula = .~. + (1|day))
#N2.month.n <- update(N2.glmm.day.plate, formula = .~. + n + factor(year), glmerControl(optimizer = "bobyqa"))
print(anova(N2.glmm.plate, N2.glmm.day.plate))
print(anova(N2.glmm.day, N2.glmm.day.plate))
print(summary(N2.glmm.day.plate))

RT <- read.csv(file.path(here::here(),"extdata","excluded_datasets","RT_rict-1.csv"))

elt2 <- filter(RT, band == "elt-2") %>%
  select(-band) %>% 
  rename(area.elt2 = area)

rict1 <-dplyr::filter(RT, band == "rict-1") %>% 
  select(-band) %>% 
  rename(area.rict1 = area)

RT <- merge(rict1,elt2) %>%
  mutate(normalized.expression = ( area.rict1 / area.elt2 ),
         sum_area = area.rict1 + area.elt2,
         expt_group = interaction(sample_group,PCR_expt))

sample_means <- RT %>% filter(genotype == "N2") %>%
  group_by(expt_group) %>% 
  summarize(mean.N2 = mean(normalized.expression)) %>% data.frame

sample_means_NIL <- RT %>% filter(genotype == "NIL") %>%
  group_by(expt_group) %>% 
  summarize(mean.N2 = mean(normalized.expression)) %>% data.frame
#%>%
normalize_RT <- RT %>%
  filter(PCR_expt == 2) %>%
  group_by(genotype,sample,sample_group,PCR_expt) %>% 
  summarise(mean_expr = mean(normalized.expression))



normalize_RT_byN2 <- RT %>%
  mutate(normalized.expression = 
           case_when(
             expt_group == 2.1 ~ normalized.expression,
             expt_group == 1.2 ~ (normalized.expression * sample_means[1,2]) / sample_means[2,2],
             expt_group == 2.2 ~ (normalized.expression * sample_means[1,2]) / sample_means[3,2] )) %>%
  group_by(genotype,sample,sample_group,PCR_expt) %>% 
  summarise(mean_expr = mean(normalized.expression))

normalize_RT_byNIL <- RT %>%
  mutate(normalized.expression = 
           case_when(
             expt_group == 2.1 ~ normalized.expression,
             expt_group == 1.2 ~ (normalized.expression * sample_means_NIL[1,2]) / sample_means_NIL[2,2],
             expt_group == 2.2 ~ (normalized.expression * sample_means_NIL[1,2]) / sample_means_NIL[3,2] )) %>%
  group_by(genotype,sample,sample_group,PCR_expt) %>% 
  summarise(mean_expr = mean(normalized.expression))



ggplot(normalize_RT, aes(x = genotype, y = mean_expr)) + 
  geom_point(alpha=0.5) +
  facet_grid(.~sample_group) + 
  coord_cartesian(ylim = c(0,4))

normalize_RT_byNIL %>% filter(PCR_expt == 2) %>% 
  ggplot(aes(x = genotype, y = mean_expr)) + 
  geom_point() +
  stat_summary(aes(y=mean_expr),fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar", width = 0.1, lwd = 0.35) +
  coord_cartesian(ylim = c(0,3))

ggplot(RT, aes(x = genotype, y = normalized.expression)) + 
  geom_point(aes(colour = factor(expt_group)))
 
ggplot(RT[RT$PCR_expt == 2,], aes(x = replicate, y = normalized.expression))  +
  geom_point(aes(colour = factor(sample))) + geom_line(aes(group = sample))

             
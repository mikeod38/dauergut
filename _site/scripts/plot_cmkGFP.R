path<-getwd()
patterns <- c("N2","ft7")
files<-lapply(patterns, function(x) {list.files(path, pattern = x)}) %>% unlist

data<-lapply(files, function(x) {data.frame(read.csv(x), name = x)}) %>% bind_rows() %>% 
  separate(name, c("genotype", "cond","cell"), sep = "_", extra = "drop") %>% rename(position = X, intensity = Y) %>%
  dplyr::filter(position < 5.01) %>%
  group_by(genotype, cond, cell) %>% 
  mutate(min = min(intensity), max = max(intensity), 
         rescaled = (intensity-min)/(max-min),
         norm.pos = ifelse(position < 2.5,position,(abs(position - 5))),
         group = interaction(genotype,cond),
         cell.id = interaction(genotype, cell)) #%>%

data$cell.id <- factor(data$cell.id)
data$genotype <- factor(data$genotype, levels = patterns)
data$group <- factor(data$group, levels = c("WT.starved", "WT.fed", "ft7.starved", "ft7.fed"))
subset <- sample(levels(data$cell.id),60)
  
data.median <- data %>% group_by(genotype,norm.pos,group,cond, cell.id) %>% summarize(median = median(rescaled))

data %>% ggplot(aes(x = norm.pos)) +
  geom_smooth(aes(y=rescaled,group = group, linetype = cond, colour = genotype), 
              method = "lm", formula = y ~ poly(x,3)) +
geom_smooth(aes(y=rescaled,group = group, linetype = cond, colour = genotype),
              method = "loess") +
  labs(
    x = "position (um)", y = "normalized mean intensity"
  ) + scale_colour_viridis(discrete = TRUE) + 
  theme_my + facet_wrap(~cell)

lm0 <- lm(data, formula = rescaled ~ genotype + cond + poly(norm.pos,3, raw = TRUE))
lm1 <- lm(data, formula = rescaled ~ genotype + cond * poly(norm.pos,3, raw = TRUE))
lm2 <- lm(data, formula = rescaled ~ cond + genotype * poly(norm.pos,3, raw = TRUE))
lm3 <- lm(data, formula = rescaled ~ genotype * cond * poly(norm.pos,3, raw = TRUE))

glmer <- glmer(data, 
               formula = rescaled ~ genotype * cond * poly(norm.pos,3, raw = TRUE) 
               + (1|cell.id))
glmer.0 <- glmer(data, 
               formula = rescaled ~ cond + genotype * poly(norm.pos,3, raw = TRUE) 
               + (1|cell.id))

stan.glmer <- stan_glmer(data, 
                    formula = rescaled ~ cond * genotype * poly(norm.pos,3, raw = TRUE) 
                    + (1|cell.id))

data.median %>% dplyr::filter(!cell.id %in% outliers) %>% ggplot(aes(x = norm.pos, group = genotype)) +
  geom_smooth(aes(y=median, group = group, linetype = cond, colour = genotype), method = "loess") + theme_solarized(light=FALSE) + 
  labs(
    x = "position (um)", y = "normalized median intensity"
  ) + scale_colour_viridis(discrete = TRUE) + theme_my

data %>% dplyr::filter(!cell.id %in% outliers) %>% group_by(genotype, cond, group, cell, cell.id) %>% 
  summarize(mean = mean(intensity), median = median(intensity)) %>%
  ggplot(aes(x = group,y=median)) + geom_boxplot()
  
  outliers <- c("ft7.2", "ft7.5L.csv","ft7.7R.csv")

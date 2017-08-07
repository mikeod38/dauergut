plot_CIs_ppt<-function (df, title, plot.contrasts, plot.contrasts.2, ypos, type, offset) {
  #generate plots for different types of data
  #all using genotype as predictor factor
  #type = "dauer" or "grid" or "GFP" or mRNA expression plot
  #set common values for all plots:
  box.width <- 0.4
  point.size <- 1
  line.width <- 0.3
  text.size <- 5
  p <- ggplot(df, aes(x=genotype)) + #x-layer
    theme_my +
    # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL),
    #               width=0,colour ="black", lwd=line.width) + #90% cred int for 95% one-sided H0
    # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
    #               width=0,colour = "darkred", lwd = line.width+0.7) + #50% cred interval
    # geom_point(data=mixed, aes(x=x.pos, y=mean), colour="black", size=point.size+0.2) + #model mean
    scale_x_discrete(labels=function(x) sub(" ","\n",x,fixed=TRUE)) +
    stat_summary(aes(x=genotype, y=ypos), geom="text", label=plot.contrasts, show.legend = TRUE, size=text.size) + # pvalues
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 15),
          axis.line = element_line(size=0.2),
          axis.title = element_text(size=16))
  if(type == "dauer") {
    p1 <- p + 
      geom_boxplot(aes(y=pct), width=box.width, outlier.shape=NA, lwd=line.width, fill="grey", notch = FALSE) +
      geom_point(aes(y=pct),size=point.size-(0.7*point.size),alpha=0.75) + 
      labs(title = title,
           y = "proportion dauer",
           x = "genotype"
      ) +
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75, 1.0))
  } else {#roaming plots
    if(type == "grid") {
      p1 <- p + geom_boxplot(aes(y=n_entries), width=box.width, outlier.shape=NA, lwd=line.width, fill = "#33CCFF",notch = FALSE) +
        #geom_point(aes(y=n_entries), size=point.size-(0.7*point.size),alpha=0.75) + 
        labs(title = title,subtitle="roaming plot",y = "grid entries",x = "genotype")
    } else {#GFP expression plots
      if(type == "GFP") {
        p1 <- p + geom_boxplot(aes(y=cell.norm),width=box.width, outlier.shape=NA, lwd=line.width, alpha = 0.75, fill="#339900", notch = FALSE) +
          #p1 <- p + geom_violin(aes(y=cell.norm), outlier.shape=NA, lwd=line.width, alpha = 0.75, fill="#339900") +
          geom_quasirandom(aes(y=cell.norm),width = 0.05, method = 'smiley', size=point.size-(0.7*point.size),alpha=0.75) + 
          labs(title = title,subtitle="GFP plot",y = "normalized expression",x = "genotype")
      } else {#mRNA plots
        p1 <- p + geom_boxplot(aes(y=cell.norm),width=box.width, outlier.shape=NA, lwd=line.width, alpha = 0.75, fill="#990000", notch = FALSE) +
          geom_quasirandom(aes(y=cell.norm),width = 0.05, method = 'smiley',size=point.size-(0.7*point.size),alpha=0.75) + 
          labs(title = title,subtitle ="mRNA FISH plot", y = "normalized expression",x = "genotype")
      }
    }
  }
  if(missing(plot.contrasts.2)) {
    p1 + stat_summary(aes(x=as.numeric(as.factor(genotype)) + 0.3, y=0),
                      fun.data = fun_length, geom = "text", size = text.size)
  } else { # add secondary comparisons
    p1 + stat_summary(aes(x=genotype, y=ypos-offset), geom="text", label=plot.contrasts.2, size=text.size, colour="red") +
      stat_summary(aes(x=as.numeric(as.factor(genotype)) + 0.3, y=0),
                   fun.data = fun_length, geom = "text", size = text.size)
  }
}


ggplot(rict1.food, aes(x=food)) +
  stat_summary(aes(y=pct, group=day), colour = "black", fun.y = mean, geom = "line", alpha = 0.2, linetype = 2) +
  geom_boxplot(width=0.3,aes(y=pct, fill=food), outlier.shape = NA, lwd=0.6, notch =FALSE) +
  geom_point(aes(y=pct), size = 0.7, alpha = 0.75) +
  labs(title = "rict-1 mutants are deficient in food suppression",
       y = "proportion dauer",
       x = "food") +
  facet_grid(.~genotype, switch="both") +
  scale_fill_manual(values = c("grey", "#FF9933")) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75, 1.0)) +
  scale_x_discrete(labels=function(x) sub(" ","\n",x,fixed=TRUE)) +
  # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL),
  #                 width=0,colour = "black", lwd=0.2) + #90% cred int for 95% one-sided H0
  # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
  #                 width=0,colour = "darkred", lwd = 0.2+0.7) + #50% cred interval
  # geom_point(data=mixed, aes(x=x.pos, y=mean), colour="black", size=1) + #model mean
  geom_text(data = plot.contrasts.H0, aes(x=2, label = prange, y = 1.075, group = NULL), size = 5) +
  theme_my_classic + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size=20))


####### daf-28 data
df %>% ggplot(aes(x=genotype, y = cell.norm)) + 
  geom_boxplot(aes(x=genotype, y=cell.norm, fill = neuron), outlier.shape = NA, lwd=0.6, width = 0.3, alpha = 0.75) +
  scale_fill_manual(values = c("#339900", "purple")) +
  geom_quasirandom(aes(y=cell.norm), width = 0.05, method = 'smiley', size=.5,alpha=0.75) +
  labs(title = "daf-28 is also reduce in rict-1 mutants",
       y = "mean GFP intensity") +
  stat_summary(aes(x=genotype, group = genotype, y = 4000), geom="text", fun.data=box_annot, label=plot.contrasts, size = 6) +
  # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL, group = genotype),
  #               width=0,colour = "black", lwd=0.6) + #90% cred int for 95% one-sided H0
  # geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75, group = genotype),
  #               width=0,colour = "darkred", lwd = 2) + #50% cred interval
  # geom_point(data=mixed, aes(x=x.pos, y=mean, group =genotype), colour="black", size=2) + #model mean
  stat_summary(aes(as.numeric(as.factor(genotype)) + 0.3,group = genotype, y = 0),
               fun.data = fun_length, geom = "text",size = 6) +
  facet_grid(.~neuron, switch="both") + theme_my_classic + 
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size=20))

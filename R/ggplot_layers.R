#' ggplot_layers
#' 
#' adds layers to ggplots
#' 
#' @section 
#' 
#' @param width optional width for median lines
#' @examples p <- ggplot(aes(x=genotype, y=pct)
#' p + add.scatter
#' @name ggplot_layers
NULL

#' @rdname ggplot_layers

add.scatter <- function() {
  geom_quasirandom(aes(y=cell.norm),colour = "#339900", cex=1,
                   width = 0.075,size=0.3,
                   method = 'smiley')}

add.median <- function(width) {
  if(missing(width)) {
    width = 0.25
  } else {
    width = width
  }
  stat_summary(aes(y=cell.norm),fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar", width = width, lwd = 0.35)
}

add.median.dauer <- function(width) {
  if(missing(width)) {
    width = 0.25
  } else {
    width = width
  }
  stat_summary(aes(y=pct),fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar", width = width, lwd = 0.35)
}

add.quartiles <- function() {
  stat_summary(aes(y=cell.norm),fun.y = median,
               fun.ymin = function(z) {quantile(z,0.25)},
               fun.ymax = function(z) {quantile(z,0.75)},
               geom = "errorbar", width = 0.15, lwd = 0.15)
}

figure.axes <- function() {
  theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 15),
        strip.text.x = ggplot2::element_blank())
}

add.n.categorical <- function() {
  stat_summary(aes(x=as.numeric(as.factor(genotype)) + 0.3, y=0),
               fun.data = fun_length, geom = "text", size = 3)
}

add.n <- function() {
  stat_summary(aes(x=temp + 0.3, y=0),
               fun.data = fun_length, geom = "text", size = 3)
}

add.Bayes.CI <- function() {
  list(geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL),
                     width=0,colour ="grey", lwd=0.15), 
       geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
                     width=0,colour = "darkgrey", lwd = 0.15+0.7),
       geom_segment(data = mixed, aes(x = x.pos-(0.009*nrow(mixed)),
                                      y = mean, xend = x.pos+(0.009*nrow(mixed)),
                                      yend = mean), colour = "darkgrey"))
}


#alt to theme classic
theme_my_classic <- ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12),legend.key = ggplot2::element_blank())

#plotting theme I use for most plots
theme_my <- ggplot2::theme_bw() + ggplot2::theme(
  axis.line        = ggplot2::element_line(colour = "black"),
  panel.grid.major = ggplot2::element_blank(), 
  panel.grid.minor = ggplot2::element_blank(),
  panel.border     = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  legend.key       = ggplot2::element_blank(), 
  axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12)
)
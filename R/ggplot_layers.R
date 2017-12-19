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

#' @export
#' @rdname ggplot_layers
#' 

add.scatter <- function() {
  geom_quasirandom(aes(y=cell.norm),colour = "#339900", cex=1,
                   width = 0.075,size=0.3,
                   method = 'smiley')}
#' @export
#' @rdname ggplot_layers

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

#' @export
#' @rdname ggplot_layers

add.mean <- function(width) {
  if(missing(width)) {
    width = 0.25
  } else {
    width = width
  }
  stat_summary(aes(y=cell.norm),fun.y = mean,
               fun.ymin = mean,
               fun.ymax = mean,
               geom = "crossbar", width = width, lwd = 0.35, colour = "white")
}

#' @export
#' @rdname ggplot_layers

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

#' @export
#' @rdname ggplot_layers

add.quartiles <- function(width) {
  if(missing(width)) {
    stat_summary(aes(y=cell.norm),fun.y = median,
                 fun.ymin = function(z) {quantile(z,0.25)},
                 fun.ymax = function(z) {quantile(z,0.75)},
                 geom = "errorbar", width = 0.15, lwd = 0.15)
  } else {
    stat_summary(aes(y=cell.norm),fun.y = median,
                 fun.ymin = function(z) {quantile(z,0.25)},
                 fun.ymax = function(z) {quantile(z,0.75)},
                 geom = "errorbar", width = width, lwd = 0.15)
  }

}

#' @export
#' @rdname ggplot_layers

figure.axes <- function() {
  theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 15),
        strip.text.x = ggplot2::element_blank())
}

#' @export
#' @rdname ggplot_layers

add.n.categorical <- function() {
  stat_summary(aes(x=as.numeric(as.factor(genotype)) + 0.3, y=0),
               fun.data = fun_length, geom = "text", size = 3)
}

#' @export
#' @rdname ggplot_layers

add.n <- function(x_axis) {
  stat_summary(aes(x=quo(x_axis) + 0.3, y=0),
               fun.data = fun_length, geom = "text", size = 3)
}

#' @export
#' @rdname ggplot_layers

add.Bayes.CI <- function() {
  list(geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL),
                     width=0,colour ="grey", lwd=0.15), 
       geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
                     width=0,colour = "darkgrey", lwd = 0.15+0.7),
       geom_segment(data = mixed, aes(x = x.pos-(0.009*nrow(mixed)),
                                      y = mean, xend = x.pos+(0.009*nrow(mixed)),
                                      yend = mean), colour = "darkgrey"))
}

#' @export
#' @rdname ggplot_layers
#alt to theme classic
theme_my_classic <- ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12),legend.key = ggplot2::element_blank())

#' @export
#' @rdname ggplot_layers
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

theme_black = function(base_size = 12, base_family = "") {
  
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_blank(),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_blank(), #(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white", face = "italic"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

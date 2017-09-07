#' box_annot
#'
#' annotate boxplots with test, based on x(factor) position
#' @param x values is the x aesthetic in ggplot2 call
#' @export
#' @examples ggplot() + stat_summary(aes(...), geom="text", fun.data=box_annot, label=list)
box_annot<-function(x) {
  return(data.frame(y=max(x) + 0.05 + 0.1*max(x)))
}

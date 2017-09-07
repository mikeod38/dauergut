#' fun_length
#'
#' annotate plots with sample size, useful for box and scatter plots
#' @param x input is the y values for which you are plotting, default is to the minimum value
#' @export
#' @examples ggplot() + stat_summary(aes(x=...) + 0.3, y=0),
#' fun.data = fun_length, geom = "text", size = 3)
fun_length<-function (x) {
  # annotate plot with sample size
  return(data.frame(y=min(x), label = paste0("(", length(x),")")))
}

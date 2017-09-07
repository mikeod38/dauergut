#' prange
#'
#' convert p values to estimates using a contrast object from tukey_ or dunnett_contrasts.R
#' @param x contrasts object, must have p.value column
#' @return an lsmeans contrast list with 'prange' column
#' @export
#' @examples contrasts <- prange(dunnett_contrasts(fit, ...))
prange <- function(x) {
  x$prange <- ifelse(x$p.value < 0.001,"***",
                     ifelse(x$p.value < 0.01, "**",
                            ifelse(x$p.value > 0.99, "p~1",
                                   ifelse(x$p.value > 0.2, sub("^","p~",round(x$p.value,digits=2)),
                                   sub("^","p~",round(x$p.value,digits=3)
                                   )))))
  return(x)
}
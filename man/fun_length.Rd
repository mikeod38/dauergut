% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_length.R
\name{fun_length}
\alias{fun_length}
\title{fun_length}
\usage{
fun_length(x)
}
\arguments{
\item{x}{input is the y values for which you are plotting, default is to the minimum value}
}
\description{
annotate plots with sample size, useful for box and scatter plots
}
\examples{
ggplot() + stat_summary(aes(x=...) + 0.3, y=0),
fun.data = fun_length, geom = "text", size = 3)
}

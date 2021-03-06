% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MM_CIs.R
\name{MM_CIs}
\alias{MM_CIs}
\title{MM_CIs}
\usage{
MM_CIs(y, predictor)
}
\arguments{
\item{y}{input if model fit ie lm, glm, glmer etc...}

\item{predictor}{optional predictor variable for comparisons, defaults to "genotype"}
}
\description{
generate Tukey-adjusted 95 pct confidence intervals for plotting on the
untransformed response scale for linear models
using lsmeans package
}
\examples{
mixed <- MM_CIs(fit, predictor = "genotype")
}

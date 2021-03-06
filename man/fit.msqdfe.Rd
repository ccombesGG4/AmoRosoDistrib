% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.msqdfe.R
\name{fit.msqdfe}
\alias{fit.msqdfe}
\title{MSQFE: Minimum Minimum squared distance Estimate based on cdf for the 4-parameter GG4 distribution}
\usage{
fit.msqdfe(y, theta)
}
\arguments{
\item{y}{the observed sample}

\item{theta}{initial values of parameters.}
}
\value{
a vector of four-parameter estimated according observation sample + info on convergence
}
\description{
Fit of univariate distributions to non-censored data by Minimum minimum squared distance estimate based on cdf.
}
\details{
See R code example in vignettes directory
}
\seealso{
\code{\link{fit.mle}},
\code{\link{fit.mkle}},
\code{\link{fit.mjse}},
\code{\link{fit.mhe}},
\code{\link{fit.msqe}},
\code{\link{fit.mwe}},
\code{\link{fit.mhdfe}},
\code{\link{fit.mwdfe}}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.mhe.R
\name{fit.mhe}
\alias{fit.mhe}
\title{MHE: Minimum Hellinger distance Estimate based on pdf for the 4-parameter GG4 distribution}
\usage{
fit.mhe(y, theta)
}
\arguments{
\item{y}{the observed sample}

\item{theta}{initial values of parameters.}
}
\value{
a vector of four-parameter estimated according observation sample + information on convergence
}
\description{
Fit of univariate distributions to non-censored data by minimum Hellinger distance estimate based on pdf.
}
\details{
See R code example in vignettes directory
}
\seealso{
\code{\link{fit.mle}},
\code{\link{fit.mkle}},
\code{\link{fit.mjse}},
\code{\link{fit.msqe}},
\code{\link{fit.mwe}},
\code{\link{fit.mhdfe}},
\code{\link{fit.mwdfe}},
\code{\link{fit.msqdfe}}.
}

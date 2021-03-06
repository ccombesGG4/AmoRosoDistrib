% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trunc.pgg4.R
\name{trucpgg4}
\alias{trucpgg4}
\title{truncated-Cumulative Distribution function of
four-parameter Generalized Gamma Distribution}
\usage{
trucpgg4(q, a, l, c, mu = 0, trt, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{q}{vector of quantile}

\item{a}{scale parameter}

\item{l}{shape parameter}

\item{c}{power shape parameter}

\item{mu}{location parameter}

\item{trt}{vector of quantile for truncated cfd}

\item{lower.tail}{Logical; if TRUE (default), probabilities are P(X≤x), otherwise, P(X>x)}

\item{log.p}{Logical; if TRUE, probabilities p are given as log(p) only if a,l and c are strictly positives}
}
\value{
probability
}
\description{
truncated-Cumulative Distribution function of
four-parameter Generalized Gamma Distribution
}

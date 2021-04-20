\name{AmoRosoDistrib-package}
\alias{AmoRosoDistrib-package}
\alias{AmoRosoDistrib}
\docType{package}

\title{
AmoRosoDistrib -- Fitting AmoRoso Family Distributions
}

\description{

\pkg{ AmoRosoDistrib } is a suite of functions for carrying out fitting from
univariate  sample (non-censoerd) to the four-parameter
Generalized Gamma Distribution (GG4) initially be described by Amoroso in 1925
who applied it to fit income rates and to model lifetime.

Versions >= 1.0-0 of this package is the first version.

For assistance using \pkg{AmoRosoDistrib} (>= 1.0-0), please see the tutorial at:

\url{https://www.jstatsoft.org/article/view/}


Primary functions available in \pkg{AmoRosoDisrib} include:

\code{gg4}: Functions for finding the density,cumulative probability distribution (cdf), quantiles and make random draws, the hazard function
and the cumulative hazard function from Amoroso Distributions called GG4.

\code{mimic}: Fit an univariate GG4 distribution functions or GEV distribution
with the minimum distance approach. Seven distance estimations are proposed
in the objective function ("Wass" for Wasserstein distance, "Chi2"
for chi-squared distance, "Burg" for Burg entropy, "NegBurg" for Negative Burg entropy,
"Euclid" for Euclidian distance, "Hell" for Hellinger distance and
"HL1N" for half l1-norm distance).

\code{logLgg4}: Log Likelihood of GG4 distribution.

\code{MaxLikEgg4}: Estimate the maximum log-likelihood of GG4 distribution
with stats::optim function.

\code{mlegg4}: Estimate the GG4 maximum log-likelihood parameters with
stats4::mle function.

\code{mean_theo, var-theo}: Functions to estimate mean and variance from
the GG4 four-parameters.

\code{moment}: Function to estimate centred moment from a univariate population
(except for the first moment which corresponds to the mean).

\code{kurtosis, skewness}: Function to compute kurtosis (tail shape)
and Skewness (measure of symmetry) of an univariate data population.


To see how to cite \pkg{AmoRosoDistrib} in publications or elsewhere,
use \code{citation ("AmoRosoDitrib")}.

}


\references{
L. Amoroso. (1925). Ricerche intorno alla curva dei redditi, \emph{Annali de Mathematica}, \bold{series 2 (1)}, 123-159.

E.W. Stacy.  (1962). A generalization of gamma distribution. \emph{Annals of Mathematical Statistics}, \bold{33}, 1187--1192.

G. E. Crooks. (2015). \emph{The amoroso distribution}, arXiv e-prints arXiv:1005.3274v2arXiv:1005.3274v2.

H. L. Harter. (1967). Maximum-likelihood estimation of the parameters of a
four- parameter generalized gamma population from complete and censored samples.
\emph{Technometrics}, \bold{9:1}, 159--165.

N. L. Johnson, S. Kotz, N. Balakrishnan. (1995).
\emph{Continuous univariate distribution}.
2cd edition John WILEY & Sons INC, New York.

E. W. Stacy. (1962). A generalization of the gamma distribution, \emph{Ann. Math.
Statist.}, \bold{33:3}, 1187--1192.

E. W. Stacy and G. A. Mihram. (1965). Parameter estimation for a generalized gamma
distribution, \emph{Technometrics}, \bold{7:3}, 349--358.

}

\keyword{ package }


#############################################################################
#   Copyright (c) 2021 Hon Keung Tony Ng and Catherine COMBES
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################
###' Parameter Initialization

#' @description Iterative method to obtain initial estimates in using the following property: If V follows a two-parameter gamma distribution with scale parameter 1 and shape parameter l, then Y = aV^(1/c)+mu follows a GGD(a,l,c,mu).
#'
#' @param data  the observed sample
#' @param lower the lower interval value of the discrete values of the power shape parameter $c$ except zero.
#'              Default value: -20
#' @param upper the upper interval value of the discrete values of the power shape parameter $c$ except zero.
#'              Default value: +20
#'@param length size of the number of discrete values of the power shape parameter $c$ except zero.
#'              Default value: 1,000
#' @param a.pos boolean pos = TRUE a>0 else a<0
#' @param p     coefficient to initialate init.mu for a>0  - Default value: 0.99
#' @param q     coefficient to initialate init.mu for a<0  - Default value: 1.01
#' @import stats
#' @return int.theta a vector of the four-parameter (a,l,c,mu) of GG4 distribution
#'
#' @details  In using a sequence of discrete values of the c parameter fixed and  the relation that X^c follows a two-parameter gamma distribution with parameters
#' a^c and l, the method of moments estimates is used to estimate a^c and l.
#' Then, for each set of estimates of theta = (a, l, c, mu), we compute the value of the likelihood,
#' and pick the set that gives the largest likelihood as initial value  four-parameter Generalized Gamma Distribution
#' So we obtain an initialisation of a vector (a,l,c,mu) as in input for the referenced method (see also paragraph)
#' See R code examples in vignettes directory.
#'
#' @seealso
#' \code{\link{fit.mle}},
#' \code{\link{fit.mkle}},
#' \code{\link{fit.mjse}},
#' \code{\link{fit.mhe}},
#' \code{\link{fit.msqe}},
#' \code{\link{fit.mwe}},
#' \code{\link{fit.mhdfe}},
#' \code{\link{fit.mwdfe}},
#' \code{\link{fit.msqdfe}}.
#'
#' @export

init.theta <- function(data, lower = -20, upper = 20, length = 1000, a.pos = TRUE, p = 0.99, q = 1.01){
  if(all(data < 0)) {
    negdata = TRUE
    data = -data
    if (a.pos)
      a.pos=F
    else a.pos = T
  }
  else  negdata = FALSE

  seqc = seq(lower, upper, length = length)
  p = p
  q = q

  distlik <- function(theta, y)
  {a <- theta[1]
  l <- theta[2]
  c <- theta[3]
  mu <- theta[4]

  fct <- -sum(log(dgg4(y, a, l, c, mu)))

  return(fct)
  }
  ###
  # in using a sequence of discrete values of $c$ fixed and
  # the relation that $X^c$ follows a two-parameter gamma distribution with parameters
  # a^c and l, the method of moments estimates is used to estimate a^c and l.
  # Then, for each set of estimates of theta = (a, l, c, mu), we compute the value of the likelihood,
  # and pick the set that gives the largest likelihood as initial value
  # four-parameter Generalized Gamma Distribution
  ###


  if (a.pos) {

    #if(min(data)<0 && p<=1){stop("When min(data)<0 and max(data)>0, \"p\" must be strictly greater than 1")}
    if(min(data)<0 && max(data)>0)  p = 1.01

    if(min(data)>0 && p>=1)
    {stop("When min(data)>0 and max(data)>0, \"p\" must be strictly smaller than 1")}
    int.mu <- min(data)*p
  }
  else{
    if(q<=1)
    {stop("\"q\" must be strictly greater than 1")}
    int.mu <- max(data)*q
  }

  xx <- (data - int.mu)

  int.theta <- numeric(4)
  largelik <- -Inf

  for (cc in seqc[!seqc == 0])
  {
    ww   <- sign(xx)*(abs(xx)^cc)

    estl <- (mean(ww)*mean(ww))/var(ww)
    esta <- sign(mean(ww))*(abs(var(ww)/mean(ww)))^(1/cc)

    theta <- c(esta, estl, cc, int.mu)

  if (estl > 0
    && esta != 0
    && !is.na(estl)
    && !is.na(esta)
    && is.finite(distlik(theta, data))
    && !is.nan(distlik(theta, data))
    && -distlik(theta, data) > largelik)
    { int.theta <- theta
      largelik <- -distlik(theta, data)
    }
  }


  if(negdata) {
    int.theta[1] <- -int.theta[1]
    int.theta[4] <- -int.theta[4]
  }
  init.res <- c(int.theta,largelik)
  return(init.res)
}


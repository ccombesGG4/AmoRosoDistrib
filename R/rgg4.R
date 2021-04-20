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
##' GG4 Random number generation
###
#
#' @description Random generation for the four-parameter Generalized Gamma Distribution
#' @param n number of observations.
#' @param a scale parameter
#' @param l shape parameter
#' @param c power shape parameter
#' @param mu location parameter
#' @param sequence if TRUE => seq for reproductive sequence else runif for alea
#' @importFrom stats runif
#' @importFrom stats qgamma
#' @return a sample generated according to a four parameter GGD
#'
#' @seealso \code{\link{dgg4}},  \code{\link{pgg4}}, \code{\link{qgg4}}, \code{\link{hgg4}}, \code{\link{chgg4}}
#'
#' #' @examples
#' ## a = 2, l = 3, c = 4 and mu = 5
#' rr = rgg4(100,2,3,4,5,sequence = TRUE)
#' c(mean(rr),sd(rr))
#' ## [1] 7.5497902 0.3718677
#'
#' @export
rgg4 <- function(n, a, l, c, mu = 0, sequence = T)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}

  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}

  scale = 1
  shape = l # shape = l with the GG4 notations

  #Y~Generalized Gamma (a,l,c,mu)
  if (sequence){
      rseq = seq(0.0,1,length = n+2)
      rseq = rseq[!rseq == 0]
      rseq = rseq[!rseq == 1]
      Y <- a*(qgamma(rseq, shape = shape, scale = scale))^(1/c) + mu
    }
    else Y <- a*(qgamma(runif(n), shape = shape, scale = scale))^(1/c) + mu
  return(Y)
}

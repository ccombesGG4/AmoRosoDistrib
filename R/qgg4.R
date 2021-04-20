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
###' GG4 Quantiles/Percentiles generation
###
#
#' @description GG4 QPG: Quantiles/Percentiles Generation of the four-parameter Generalized Gamma Distribution
#' @param p vector or quantiles
#' @param a scale parameter
#' @param l shape parameter
#' @param c power shape parameter
#' @param mu location parameter
#' @param lower.tail indicate lower tail. The default value is set to TRUE
#' @importFrom stats qgamma
#' @return a vector of quantile generated according to a four parameter GGD
#'
#' @seealso \code{\link{dgg4}}, \code{\link{pgg4}}, \code{\link{rgg4}},\code{\link{hgg4}}, \code{\link{chgg4}}
#'
#' @examples
#' set.seed(963)
#' qq = qgg4(runif(100),2,3,4,5)
#' c(mean(qq),sd(qq))
#' ## [1] 7.505230 0.368037
#'
#' # ----
#' # Test
#' # ----
#'
#' # a = 1, l = 2, c = -0.8 and mu = 0 (default value)
#' p <- (1:9)/10
#' pgg4(qgg4(p, 1, 2, -0.8), 1, 2, -0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#'
#' @export
qgg4 <- function(p, a, l, c, mu = 0, lower.tail = TRUE)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}

  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}

  #Y~ GenGamma(a,l,c,mu)
  Y <- a*qgamma(p,scale = 1,shape = l)^(1/c) + mu
  if (!lower.tail) Y <- 1 - Y
  return(Y)
}







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
###' GG4 Cumulative Distribution Function (cdf)
###
#' @description GG4 CDF: Cumulative distribution function the four-parameter Generalized Gamma Distribution
#' @param q vector of quantile
#' @param a scale parameter
#' @param l shape parameter
#' @param c power shape parameter
#' @param mu location parameter
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X≤x), otherwise, P(X>x)
#' @param log.p Logical; if TRUE, probabilities p are given as log(p) only if a,l and c are strictly positives
#' @import stats
#' @importFrom stats pgamma
#' @return probability
#'
#' @seealso \code{\link{dgg4}}, \code{\link{qgg4}}, \code{\link{rgg4}},\code{\link{hgg4}}, \code{\link{chgg4}}
#'
#' @examples
#' # a = 1, l = 2, c = -0.8 and mu = 0 (default value)
#' p <- (1:9)/10
#' pgg4(qgg4(p, 1, 2, -0.8), 1, 2, -0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#'
#' @export



pgg4 <- function(q, a, l, c, mu = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}

  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}


  #if (!log.p) {

     W <- ((q - mu)/a)^c
     Y = pgamma(sort(W), shape = l, lower.tail = lower.tail, log.p = log.p)

 if (!lower.tail) Y <- 1 - Y

  return(Y)
}


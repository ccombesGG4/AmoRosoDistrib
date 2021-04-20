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
###' truncated-Cumulative Distribution function of
###' four-parameter Generalized Gamma Distribution
###
#' @param q vector of quantile
#' @param trt vector of quantile for truncated cfd
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
#' @export



trucpgg4 <- function(q, a, l, c, mu = 0, trt, lower.tail = TRUE, log.p = FALSE)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}
  if (missing(mu)) mu = 0
  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}

  if (lower.tail) {
    trp <- pgg4(sort(trt), a, l, c, mu, lower.tail = TRUE, log.p = FALSE)
    Y   <- (pgg4(sort(q), a, l, c, mu, lower.tail = TRUE, log.p = FALSE) - trp )/(1 - trp)
    if (log.p) Y <- log(Y)

  } else {
    Y   <- pgg4(sort(q), a, l, c, mu, lower.tail= FALSE, log.p = TRUE) - pgg4(sort(trt), a, l, c, mu, lower.tail = FALSE, log.p = TRUE)
    if (!log.p) Y <- exp(Y)
  }


  return(Y)
}


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
###' GG4 Probability Density Function
###
#
#' @description GG4 PDF: Probability Density Function of the four-parameter Generalized Gamma Distribution
#' @param x vector of quantile
#' @param a scale parameter
#' @param l shape parameter
#' @param c power shape parameter
#' @param mu location parameter
#' @param log.p Logical; if TRUE, probabilities p are given as log(p)
#' @import stats
#' @return density
#'
#' @seealso \code{\link{pgg4}}, \code{\link{qgg4}}, \code{\link{rgg4}},\code{\link{hgg4}}, \code{\link{chgg4}}
#'
#' @examples
#' set.seed(963)
#' rr = rgg4(100,2,3,4,5,sequence = FALSE)
#' c(mean(rr),sd(rr))
#' ## [1] 7.505230 0.368037
#' q = sort(rr)
#' dd = dgg4(q,2,3,4,5)
#' plot(dd)
#'
#' @export

dgg4 <- function(x,a,l,c,mu = 0, log.p=F)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}

  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}

  # 1/gamma(l)* |c/a| * [(y-mu)/a]^(lc-1) * exp[-((y-mu)/a)^c ]
  Y <- 1/gamma(l) * abs(c/a) *((x - mu)/a)^(l*c - 1)*exp(-((x - mu)/a)^c)

  if (log.p) Y = log(Y)
  return(Y)
}

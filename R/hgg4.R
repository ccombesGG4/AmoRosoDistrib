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
###' GG4 Hazard Function
###
#
#' @description GG4 HF: Hazard function of the four-parameter Generalized Gamma Distribution
#' @param x data
#' @param a scale parameter
#' @param l shape parameter
#' @param c power shape parameter
#' @param mu location parameter
#' @return hazard ratio
#'
#' @seealso \code{\link{dgg4}},  \code{\link{pgg4}}, \code{\link{qgg4}}, \code{\link{rgg4}}, \code{\link{chgg4}}
#'
#' @examples
#' ## a = 2, l = 3, c = 4 and mu = 5
#' rr = rgg4(100,2,3,4,5,sequence = TRUE)
#' c(mean(rr),sd(rr))
#' ## [1] 7.5497902 0.3718677
#' q = sort(rr)
#' hh = hgg4(q,2,3,4,5)
#' summary(hh)
#' ## Min.    1st Qu.  Median    Mean    3rd Qu.  Max.
#' ## 0.07809 1.15079  2.06266   2.27508 3.16877  6.60587
#'
#' @export
hgg4 <- function(x,a,l,c,mu = 0)
{
  if (missing(l)) {stop("shape parameter \"l\" not given")}
  if (missing(c)) {stop("power shape parameter \"c\" not given")}
  if (missing(a)) {stop("scale parameter \"a\" not given")}

  if (l < 0)  {stop("Shape parameter \"l\" must be strictly positive")}
  if (a == 0) {stop("Scale parameter \"a\" should be different from zero")}
  dgg4(x = x, a = a, l = l, c = c, mu = mu) /
    pgg4(q = x, a = a, l = l, c = c, mu = mu, lower.tail = FALSE)
}

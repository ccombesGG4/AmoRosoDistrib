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
###' upgamma function
###'
###
#' @param x vector of quantile
#' @param shape shape parameter
#' @import stats
#' @importFrom stats pgamma
#' @export

upgamma <- function(shape, x) {
  if (missing(shape)) {stop("shape parameter \"shape\" not given")}
  if (shape < 0)  {stop("Shape parameter \"shape\" must be strictly positive")}
  pgamma(x, shape=shape, lower.tail = FALSE)*gamma(shape)
}


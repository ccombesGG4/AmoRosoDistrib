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
##' Skewness of a data population (measure of symmetry)
###
###'
###
#

#' @param data a sample (set of observations)
#' @param type a character string coding -- "sample" is computed from mean and standard deviation
#' of sample. "fisher" corresponds to fisher'skewness formula.
#' @return S the skewness from a sample
#'
#' @export

skewness <- function(data, type = "sample"){
  if (!is.element(type, c("sample", "fisher")))
  stop(paste("The ", type, "does not exist. Character string coding function must be
             sample ou fisher"))
  n = length(data)
  if (type == "sample")
    S = (n/((n - 1)*(n - 2)))*sum(((data - mean(data))/sqrt(var(data)))^3)
  if (type == "fisher")
    S = (((n*(n - 1)^0.5)/(n - 2))*sum(data^3/n))/(sum(data^2/n)^3/2)
  return(S)
}

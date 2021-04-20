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
###' To compute Bias-reduced estimate of mu
#' @param est_theta list of the 4-parameter GG4 with Est_theta[1]= a, Est_theta[2]=l, Est_theta[3]=C,and Est_theta[4]=mu
#' @param y data corresponding to the observed population

#' @import stats
#' @return bre_mu  Bias-reduced estimate of mu
#'
#' @export

mu_bias_reduced_estimate <- function(est_theta,y)
{ a <- est_theta[1]
  l <- est_theta[2]
  c <- est_theta[3]

  if(!is.na(a)){

## Estimated bias of mu:
  estbias <- est_theta[1]*gamma(est_theta[2] + 1/est_theta[3])/gamma(est_theta[2])

 ## Bias-reduced estimate of mu
  if (a>0)
    bre_mu <- min(min(y), mean(y) - estbias)
  else
    bre_mu <- max(max(y), max(y) - estbias)

  return(bre_mu)
  }
  else return("NA")

}

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
##' Moment estimated from the Generalized Gamma four-parameters. the function
##' allows us to verify the results regarding parameters estimated with mimic function
###
#

#' @param theta vector of 4-parameters.
#' @param k order of the moment to be computed from the GGD parameter
#' @return moment from parameters
#' @export

  moment_theo <- function(theta,k){
    a <- theta[1]
    l <- theta[2]
    c <- theta[3]
    mu <- theta[4]

    if (k == 1) {
      E = (a*gamma(l + 1/c)) / gamma(l) + mu
      }
    else{
      #if ((abs(c)*l) > k){
        E = (a^k*gamma((c*l + k)/c)) / (gamma(l))
    }
    return(E)
 }


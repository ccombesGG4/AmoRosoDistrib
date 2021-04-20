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
##' Variance estimated from the Generalized Gamma four-parameters. the function
##' allows us to verify the results regarding parameters estimated with mimic function
###
#

#' @param theta vector of 4-parameters.
#' @return variance from parameters
#' @export

  var_theo <- function(theta){
    a <- theta[1]
    l <- theta[2]
    c <- theta[3]
    mu <- theta[4]

    V = (a^2*(gamma((c*l + 2)/c)*gamma(l) - (gamma((c*l + 1)/c))^2)) / (gamma(l))^2

    return(V)
 }


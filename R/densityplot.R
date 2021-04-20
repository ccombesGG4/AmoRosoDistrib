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
#########################################################################################
#'  Density plots observations versus model GG4 estimated
#'
##
#
#' @param y the observed sample
#' @param theta Vector of salues of the four-parameters.
#' @return a QQ-plot, densities and hystogramm of observations vs estimated model of GG4
#' @import graphics
#' @import transport
#' @import tidyverse
#' @import grDevices
#' @import ggplot2
#' @export


densityplot <- function(y, theta){

   par(mfrow = c(1,1))
   hist(y, prob = T, main = "Densities observed data vs generated data set from estimated parameters")
   lines(density(y), lwd = 4, col = "red")
   xx <- seq(min(y), max(y), length = 10000)
   par(mfrow = c(1,1))
   points(xx, dgg4(xx, theta[1], theta[2], theta[3], theta[4]), type = "l", lwd = 2, col = "blue", lty = 2)

   legend("topleft", legend = c("Obs. data", "GG4 Model"), lwd = 2, col = c("red","blue"), lty = 1:2)
}


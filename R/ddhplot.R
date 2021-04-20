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
#' Densities and histogramm observations versus model GG4 estimated
#'
##
#
#' @param y the observed sample
#' @param theta alues of parameters.
#' @param seed seed initialization
#' @return a QQ-plot, densities and hystogramm of observations vs estimated model of GG4
#' @importFrom stats density
#' @importFrom stats optim
#' @importFrom grDevices rgb
#' @import graphics
#' @import transport
#' @import tidyverse
#' @import ggplot2
#' @import UsingR
#' @export


ddhplot <- function(y, theta, seed = 986){
   #require(ggplot2)
   obs = sort(y)
   set.seed(seed)
   yest = rgg4(length(y), theta[1],theta[2],theta[3],theta[4])
   theo = sort(yest)
   min = min(min(obs),min(theo))
   max = max(max(obs),max(theo))
   par(mfrow = c(1,2))
   UsingR::QQplot(obs,theo)
   abline(0,1)
   #Histogram
   hist(y, xlim=c(min,max), col="red")
   hist(yest, add=T, col=rgb(0, 1, 0, 0.5))
}


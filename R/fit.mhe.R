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
##############################################################################################
#' MHE: Minimum Hellinger distance Estimate based on pdf for the 4-parameter GG4 distribution
#'
##
#
#' @description Fit of univariate distributions to non-censored data by minimum Hellinger distance estimate based on pdf.
#' @param y the observed sample
#' @param theta initial values of parameters.

#' @return a vector of four-parameter estimated according observation sample + information on convergence
#' @importFrom stats density
#' @importFrom stats optim
#'
#' @details See R code example in vignettes directory
#'
#' @seealso
#' \code{\link{fit.mle}},
#' \code{\link{fit.mkle}},
#' \code{\link{fit.mjse}},
#' \code{\link{fit.msqe}},
#' \code{\link{fit.mwe}},
#' \code{\link{fit.mhdfe}},
#' \code{\link{fit.mwdfe}},
#' \code{\link{fit.msqdfe}}.
#'
#' @export

fit.mhe <- function(y, theta){
  if(all(y<0)) {
    negdata = TRUE
    theta[1] = -theta[1]
    theta[4] = -theta[4]
    y = -y}
  else negdata = FALSE

  disth <- function(theta, y){
          a <- theta[1]
          l <- theta[2]
          c <- theta[3]
          mu <- theta[4]

          if (a>0) {
            if(mu<0)
              shift = mu*0.99
            else shift = mu*1.01
            dfit <- density(y, from = shift)
            }

          else {
              shift = mu*0.99
              dfit <- density(y, to = shift)
            }


          xx <- dgg4(dfit$x, a, l, c, mu)
          yy <- dfit$y

          fct = sum((sqrt(xx) - sqrt(yy))^2)

          return(fct)
  }
  if (theta[1] < 0){
    limit = max(y)
    ui = matrix(c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), 3, 4, byrow = T)
  } else {
    limit = -min(y)
    ui = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), 3, 4, byrow = T)
  }



  mhdest <- tryCatch(constrOptim(theta = theta, f = disth, y = y,
                                 grad=NULL, method="Nelder-Mead",  ui=ui,
                                 ci = c(0, 0, limit), outer.iterations=1000,outer.eps=1e-10), error=function(e) {NA})
      if (negdata) {
        mhdest[1]$par =  -mhdest[1]$par
        mhdest[4]$par =  -mhdest[4]$par
      }



  return(mhdest)
   }


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

#'  MWDFE: Minimum Wasserstein distance Estimate based on cdf for the 4-parameter GG4 distribution
#'
##
#' @description Fit of univariate distributions to non-censored data by minimum Wasserstein distance estimate based on cdf.
#' @param y the observed sample
#' @param theta initial values of parameters.
#' @param d exponent of the Wasserstein metric
#' @return a vector of four-parameter estimated according observation sample + info on convergence
#' @importFrom stats density
#' @importFrom stats optim
#' @import transport
#'
#' @seealso
#' \code{\link{fit.mle}},
#' \code{\link{fit.mkle}},
#' \code{\link{fit.mjse}},
#' \code{\link{fit.mhe}},
#' \code{\link{fit.msqe}},
#' \code{\link{fit.mwe}},
#' \code{\link{fit.mhdfe}},
#' \code{\link{fit.msqdfe}}.
#'
#' @details See R code example in vignettes directory
#'
#' @export


fit.mwdfe <- function(y, theta, d = 1){
        distwddf <- function(theta, y)  {
          a  <- theta[1]
          l  <- theta[2]
          c  <- theta[3]
          mu <- theta[4]

          npfit <- (1:length(y) - 0.5)/length(y)
          #fct <- transport::wasserstein1d(pgg4(sort(y), a, l, c, mu), npfit, d)
          fct <- mean(abs(pgg4(sort(y), a, l, c, mu) - npfit)^d)^(1/d)
          return(fct)
        }

        if (theta[1] < 0){
                limit = max(y)
                ui = matrix(c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), 3, 4, byrow = T)
        } else {
                limit = -min(y)
                ui = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), 3, 4, byrow = T)
        }


        mswdfest <- constrOptim(theta = theta, f = distwddf, y = y,
                                 grad=NULL, method="Nelder-Mead",  ui=ui,
                                 ci = c(0, 0, limit), outer.iterations=1000,outer.eps=1e-10)


       return(mswdfest)
   }


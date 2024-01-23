### sCorrect-getVarCov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:00) 
## Version: 
## Last-Updated: jan 23 2024 (10:25) 
##           By: Brice Ozenne
##     Update #: 85
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov2
#' @title Residual Variance-Covariance Matrix With Small Sample Correction.
#' @description Reconstruct the residual variance-covariance matrix from a latent variable model. 
#' It is similar to \code{nlme::getVarCov} but with small sample correction.
#' @name getVarCov2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#' 
#' @return A matrix with as many rows and column as the number of endogenous variables
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the residuals.
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' n <- 101
#'
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' Id <- findInterval(runif(n), seq(0.1,1,0.1))
#' data.df <- rbind(data.frame(Y=Y1,G="1",Id = Id),
#'            data.frame(Y=Y2,G="2",Id = Id)
#'            )
#'
#' #### latent variable models ####
#' library(lava)
#' e.lvm <- estimate(lvm(Y ~ G), data = data.df)
#' getVarCov2(e.lvm)
#' 
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`getVarCov2` <-
    function(object, ...) UseMethod("getVarCov2")

## * getVarCov2.lvmfit
#' @rdname getVarCov2
#' @export
getVarCov2.lvmfit <- function(object, ssc = lava.options()$ssc, ...){

    return(getVarCov2(estimate2(object, ssc = ssc, ...)))

}

## * getVarCov2.lvmfit2
#' @rdname getVarCov2
#' @export
getVarCov2.lvmfit2 <- function(object, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }
    

    Omega <- object$sCorrect$moment$Omega
    attr(Omega, "detail") <- NULL
    return(Omega)
    
}

######################################################################
### sCorrect-getVarCov2.R ends here

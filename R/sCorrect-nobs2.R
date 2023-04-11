### sCorrect-nobs2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan  4 2022 (14:37) 
## Version: 
## Last-Updated: Jan  6 2022 (15:53) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - nobs2
#' @title Effective Sample Size.
#' @description Extract the effective sample size, i.e. sample size minus the loss in degrees of freedom caused by the estimation of the parameters.
#' @name nobs2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object.
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the leverage.
#' 
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#' 
#' @return Numeric vector of length the number of endogenous variables.
#' 
#' @concept extractor
#' @keywords smallSampleCorrection
#' 
#' @export
`nobs2` <-
    function(object, ssc, ...) UseMethod("nobs2")

## * nobs2.lvmfit
#' @rdname nobs2
#' @export
nobs2.lvmfit <- function(object, ssc = lava.options()$ssc, ...){

    return(nobs2(estimate2(object, ssc = ssc, ...)))

}

## * nobs2.lvmfit2
#' @rdname nobs2
#' @export
nobs2.lvmfit2 <- function(object, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    return(stats::nobs(object) - colSums(leverage2(object, format = "wide")))
}

##----------------------------------------------------------------------
### sCorrect-nobs2.R ends here

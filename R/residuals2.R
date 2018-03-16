### residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:05) 
## Version: 
## Last-Updated: mar 12 2018 (17:51) 
##           By: Brice Ozenne
##     Update #: 912
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - residuals2
#' @title Extract Corrected Residuals
#' @description Extract correct residuals from a gaussian linear model.
#' @name residuals2
#' 
#' @param object a \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} object.
#' @param param [optional] the fitted parameters.
#' @param data [optional] the data set.
#' @param ... arguments to be passed to \code{sCorrect}.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the residuals.
#' 
#' @return a matrix containing the residuals relative to each sample (in rows)
#' and each endogenous variable (in column).
#'
#' @examples
#' ## simulate data
#' set.seed(10)
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#' d <- sim(m,20, latent = FALSE)
#'
#' ## standard linear model
#' e.lm <- lm(Y1~Y2, data = d)
#' sCorrect(e.lm) <- TRUE
#' 
#' sigma(e.lm)^2
#' mean(residuals(e.lm)^2)
#' mean(residuals2(e.lm)^2)
#' 
#' ## latent variable model
#' e.lvm <- estimate(m, data = d)
#' sCorrect(e.lvm) <- TRUE
#' mean(residuals2(e.lvm)^2)
#'
#' @concept small sample inference
#' @export
`residuals2` <-
    function(object, ...) UseMethod("residuals2")

## * residuals2.lm
#' @rdname residuals2
#' @export
residuals2.lm2 <- function(object, param = NULL, data = NULL, ...){

    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        args$score <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }
    return(object$sCorrect$epsilon)   
}

## * residuals2.gls
#' @rdname residuals2
#' @export
residuals2.gls2 <- residuals2.lm2

## * residuals2.lme
#' @rdname residuals2
#' @export
residuals2.lme2 <- residuals2.lm2

## * residuals2.lvmfit
#' @rdname residuals2
#' @export
residuals2.lvmfit2 <- residuals2.lm2


##----------------------------------------------------------------------
### residuals2.R ends here


### lTest.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: feb  5 2018 (10:49) 
##           By: Brice Ozenne
##     Update #: 669
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - lTest
#' @title Test linear hypotheses using Wald statistic
#' @description Test linear hypotheses using Wald statistic. Deprecated see compare2.
#' @name lTest
#'
#' @param object a \code{lvm} object.
#' @param C [optional] a contrast matrix.
#' @param ... arguments to be passed to \code{compare2}.
#' @export
`lTest` <-
  function(object, ...) UseMethod("lTest")

## * lTest.lm
#' @rdname lTest
#' @export
lTest.lm <- function(object, C = NULL, ...) {
    .Deprecated("compare2")
    if(!is.null(C)){
        warning("argument C is deprecated, use contrast instead \n")
        
    }
    return(compare2(object = object, ...))
}

## * lTest.gls
#' @rdname lTest
#' @export
lTest.gls <- lTest.lm

## * lTest.lme
#' @rdname lTest
#' @export
lTest.lme <- lTest.lm

## * lTest.lvmfit
#' @rdname lTest
#' @export
lTest.lvmfit <- lTest.lm





##----------------------------------------------------------------------
### lTest.R ends here

### sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan 22 2022 (13:50) 
## Version: 
## Last-Updated: Jan 22 2022 (14:05) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Depreciated Method For Small Sample Correction
#' @description Depreciated method for small sample correction, now replaced by the \code{\link{estimate2}} method.
#' @name sCorrect
#' 
#' @param object,x a \code{lvmfit} object.
#' @param value not used.
#' @param ... not used.
#'

#' @rdname sCorrect
#' @export
`sCorrect` <-
    function(object, ...) UseMethod("sCorrect")

#' @rdname sCorrect
#' @export
`sCorrect.default` <- function(object, ...){
    .Defunct("estimate2", package = "lavaSearch2", msg = "Function sCorrect has been removed from lavaSearch2 version 2.0.0 or higher and replace by the estimate2 method. \n")
}

#' @rdname sCorrect
#' @export
`sCorrect<-` <-
  function(x, ..., value) UseMethod("sCorrect<-")

#' @rdname sCorrect
#' @export
`sCorrect<-.default` <- function(x, ..., value){
    .Defunct("estimate2", package = "lavaSearch2", msg = "Function sCorrect has been removed from lavaSearch2 version 2.0.0 or higher and replace by the estimate2 method. \n")
}

##----------------------------------------------------------------------
### sCorrect.R ends here

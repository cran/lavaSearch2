### sCorrect-getIndexOmega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 25 2019 (10:52) 
## Version: 
## Last-Updated: jan 17 2022 (14:26) 
##           By: Brice Ozenne
##     Update #: 103
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Identify the Endogenous Variables
#' @description Identify the endogenous variables, i.e., returns a vector with length the number of observations,
#' whose values are the index of the repetitions.
#' @name getIndexOmega
#'
#' @param object a \code{lvmfit} object.
#' @param data dataset.
#' @param ... [internal] Only used by the generic method.
#'  
#' @concept extractor
#' @keywords internal
`.getIndexOmega` <-
    function(object, data, ...) UseMethod(".getIndexOmega")

## * Examples
#' @examples
#' \dontrun{
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### lvm model ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' ## lavaSearch2:::.getIndexOmega(e.lvm, data = dW)
#' }

## * .getIndexOmega.lvm
#' @rdname getIndexOmega
.getIndexOmega.lvm <- function(object, data, ...){

    ## ** check missing value in exogenous variables
    name.exogenous <- exogenous(object)
    missing.var <- name.exogenous[name.exogenous %in% names(data) == FALSE]

    if(length(missing.var)>0){
        cat2bin <- var2dummy(object$model, var = names(data), data = data)
        name.exogenous[name.exogenous %in% missing.var] <- names(cat2bin)[cat2bin %in% missing.var]
        name.exogenous <- unique(name.exogenous)
    }
    test.na <- Reduce("+",lapply(name.exogenous, function(iCol){is.na(data[[iCol]])}))
    if(any(test.na>0)){
        stop("Does not support missing values in exogenous variables. \n",
             "Consider removing the corresponding rows in the dataset. \n")
    }

    ## ** index.Omega
    n.obs <- NROW(data)
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)

    M.index <- matrix(1:n.endogenous, nrow = n.obs, ncol = n.endogenous, byrow = TRUE)
    index.na <- which(is.na(as.data.frame(data)[,name.endogenous]))
    if(length(index.na)>0){
        M.index[index.na] <- NA
    }
    return(as.double(t(M.index)))
}

## * .getIndexOmega.lvmfit
#' @rdname getIndexOmega
.getIndexOmega.lvmfit <- .getIndexOmega.lvm

######################################################################
### sCorrect-getIndexOmega.R ends here

### sampleRepeated.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (15:37) 
## Version: 
## Last-Updated: Jan 11 2022 (17:39) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Simulate Repeated Measurements over time
#' @description Simulate repeated measurements over time (one factor model).
#' @name sampleRepeated
#' 
#' @param n [integer] sample size.
#' @param n.Xcont [integer] number of continuous covariates acting on the latent variable.
#' @param n.Xcat [integer] number of categorical covariates acting on the latent variable.
#' @param n.rep [integer] number of measurement of the response variable.
#' @param format [character] should the dataset be returned in the \code{"long"} format or in the \code{"wide"} format.
#'
#' @return a \code{data.frame} object.
#'
#' @examples
#' 
#' sampleRepeated(10, format = "wide")
#' sampleRepeated(10, format = "long")

## * sampleRepeated
#' @rdname sampleRepeated
#' @export
sampleRepeated <- function(n, n.Xcont = 2, n.Xcat = 2, n.rep = 5, format = "long"){

    ## ** check arguments
    format <- match.arg(format, choices = c("wide","long"))
    
    ## ** define lvm
    m <- lava::lvm()
    idvars <- "id"
    distribution(m, ~id) <- function(n, ...){return(1:n)}    
    for(iY in 1:n.rep){ ## iY <- 1
        regression(m) <- stats::as.formula(paste0("Y",iY,"~eta"))
    }
    if(n.Xcont>0){
        for(iXcont in 1:n.Xcont){ ## iY <- 1
            regression(m) <- stats::as.formula(paste0("eta~X",iXcont))
        }
        idvars <- c(idvars,paste0("X",1:n.Xcont))
    }
    
    if(n.Xcat>0){
        for(iXcat in 1:n.Xcat){ ## iY <- 1
            regression(m) <- stats::as.formula(paste0("eta~Z",iXcat))
            categorical(m, labels=c("a","b","c")) <- paste0("Z",iXcat)
        }
        idvars <- c(idvars,paste0("Z",1:n.Xcat))
    }    
    latent(m) <- ~eta

    dW <- lava::sim(m, n = n, latent = FALSE)
    if(format == "wide"){
        return(dW)
    }else{
        dL <- stats::reshape(dW,
                             direction = "long",
                             idvar = idvars,
                             varying = list(paste0("Y",1:n.rep)),
                             v.names = "Y",
                             timevar = "time")
        rownames(dL) <- NULL
        return(dL)
    }
    
}

######################################################################
### sampleRepeated.R ends here

### sCorrect-information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (14:17) 
## Version: 
## Last-Updated: Jan 17 2022 (22:54) 
##           By: Brice Ozenne
##     Update #: 445
##----------------------------------------------------------------------
## 
### Commentary: 
## Compute information, hessian, and first derivative of information
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - information2
#' @title  Expected Information With Small Sample Correction. 
#' @description  Extract the expected information matrix from a latent variable model.
#' Similar to \code{lava::information} but with small sample correction.
#' @name information2
#'
#' @param object,x a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the information matrix.
#'
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#'
#' @return A matrix with as many rows and columns as the number of coefficients.
#' 
#' @examples
#' #### simulate data ####
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- Sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#' 
#' #### linear models ####
#' e.lm <- lm(formula.lvm,data=d)
#' 
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' information(e.lvm)
#' information2(e.lvm)
#' information2(e.lvm)[1:4,1:4] -  solve(vcov(e.lm))
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`information2` <-
  function(object, as.lava, ssc, ...) UseMethod("information2")

## * information2.lvmfit
#' @rdname information2
#' @export
information2.lvmfit <- function(object, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(information(estimate2(object, ssc = ssc, ...), as.lava = as.lava))

}

## * information2.lvmfit2
#' @rdname information2
#' @export
information2.lvmfit2 <- function(object, as.lava = TRUE, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    out <- object$sCorrect$information[names(object$sCorrect$skeleton$originalLink2param),
                                       names(object$sCorrect$skeleton$originalLink2param),
                                       drop=FALSE]
    if(as.lava==FALSE){
        dimnames(out) <- list(as.character(object$sCorrect$skeleton$originalLink2param),
                              as.character(object$sCorrect$skeleton$originalLink2param))
    }
    return(out)
}

## * information.lvmfit2
#' @rdname information2
#' @export
information.lvmfit2 <- function(x, ...){ ## necessary as first argument of information must be x 
    information2(x, ...)
}

## * .information2
#' @title Compute the Expected Information Matrix From the Conditional Moments
#' @description Compute the expected information matrix from the conditional moments.
#' @name information2-internal
#' 
#' @details \code{calc_information} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.information2 <- function(dmu, dOmega, OmegaM1,
                          missing.pattern, unique.pattern, name.pattern,
                          grid.mean, grid.var, name.param,
                          leverage, weights = NULL, n.cluster){
    if(lava.options()$debug){cat(".information2\n")}
    if(is.null(weights)){weights <- rep(1,n.cluster)}
    
    ## ** Prepare
    n.grid.mean <- NROW(grid.mean)
    n.grid.var <- NROW(grid.var)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    Info <- matrix(0, nrow = n.param, ncol = n.param,
                   dimnames = list(name.param,name.param))
    if(length(dmu)>0){
        index.mean <- 1:n.grid.mean
    }else{
        index.mean <- NULL
    }
    if(length(dOmega)>0){
        index.var <- 1:n.grid.var
    }else{
        index.var <- NULL
    } 

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        if(!is.null(leverage)){
            iN.corrected <- sum(weights[iIndex]) - colSums(leverage[iIndex,iY,drop=FALSE])
        }else{
            iN.corrected <- sum(weights[iIndex]) 
        }
        
        ## *** Information relative to the mean parameters
        for(iG in index.mean){ # iG <- 1
            iP1 <- grid.mean[iG,1]
            iP2 <- grid.mean[iG,2]
            Info[iP1,iP2] <- Info[iP1,iP2] + sum(rowSums(dmu[[iP1]][iIndex,iY,drop=FALSE] %*% iOmegaM1 * dmu[[iP2]][iIndex,iY,drop=FALSE])*weights[iIndex])
        }

        ## *** Information realtive to the variance parameters
        for(iG in index.var){ # iG <- 2
            iP1 <- grid.var[iG,1]
            iP2 <- grid.var[iG,2]

            ## NOTE: normally tr(ABAC)=tr(ACAB) but because of the factor n.correct this is no more the case
            ##       so the information may slightly dependent on the ordering of the parameters
            iDiag <- diag(iOmegaM1 %*% dOmega[[iP1]][iY,iY,drop=FALSE] %*% iOmegaM1 %*% dOmega[[iP2]][iY,iY,drop=FALSE])
                        
            Info[iP1,iP2] <- Info[iP1,iP2] + 1/2*sum(iDiag*iN.corrected)
        }        
    }

    ## ** Make Info a symmetric matrix
    Info <- symmetrize(Info, update.upper = NULL)
        
    ## ** export
    return(Info)
}






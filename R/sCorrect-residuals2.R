### sCorrect-residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:17) 
## Version: 
## Last-Updated: jan 23 2024 (10:26) 
##           By: Brice Ozenne
##     Update #: 141
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Residuals With Small Sample Correction.
#' @description Extract residuals from a latent variable model.
#' Similar to \code{stats::residuals} but with small sample correction.
#' @name residuals2
#' 
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param type [character] the type of residual to extract:
#' \code{"response"} for raw residuals,
#' \code{"studentized"} for studentized residuals,
#' \code{"normalized"} for normalized residuals.
#' @param format [character] Use \code{"wide"} to return the residuals in the wide format (one row relative to each sample).
#' Otherwise use \code{"long"} to return the residuals in the long format.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the residuals.
#'
#' The raw residuals are defined by  observation minus the fitted value:
#' \deqn{
#' \varepsilon = (Y_1 - \mu_1, ..., Y_m - \mu_m)
#' }
#' The studentized residuals divided the raw residuals relative to each endogenous variable by the modeled variance of the endogenous variable.
#' \deqn{
#' \varepsilon_{stud} =(\frac{Y_1 - \mu_1}{\sigma_1}, ..., \frac{Y_m - \mu_m}{\sigma_m})
#' }
#' The normalized residuals multiply the raw residuals by the inverse of the square root of the modeled residual variance covariance matrix.
#' \deqn{
#' \varepsilon_{norm} = \varepsilon \Omega^{-1/2}
#' }
#' @return a matrix containing the residuals relative to each sample (in rows)
#' and each endogenous variable (in column).
#' 
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`residuals2` <-
    function(object, type, format, ssc, ...) UseMethod("residuals2")

## * Examples
#' @rdname residuals2
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
#' residuals(e.lvm)
#' residuals2(e.lvm)
#' residuals(e.lvm) - residuals2(e.lvm)
#'

## * residuals2.lvmfit
#' @export
residuals2.lvmfit <- function(object, type = "response", format = "wide", ssc = lava.options()$ssc, ...){

    return(residuals(estimate2(object, ssc = ssc, ...), type = type, format = format))

}

## * residuals2.lvmfit2
#' @export
residuals2.lvmfit2 <- function(object, type = "response", format = "wide", ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    format <- match.arg(format, choices = c("long","wide"))

    residuals <- .normalizeResiduals(epsilon = object$sCorrect$residuals,
                                     Omega = object$sCorrect$moment$Omega,
                                     type = type,
                                     missing.pattern = object$sCorrect$missing$pattern,
                                     unique.pattern = object$sCorrect$missing$unique.pattern,
                                     Omega.missing.pattern = object$sCorrect$moment$Omega.missing.pattern)
    if(format == "wide"){
        return(residuals)
    }else if(format == "long"){
        endogenous <- colnames(residuals)
        n.endogenous <- length(endogenous)
        
        residualsW <- data.frame(1:NROW(residuals), residuals)
        names(residualsW) <- c("cluster",endogenous)
        residualsL <- stats::na.omit(stats::reshape(residualsW,
                                                    idvar = "cluster",
                                                    direction = "long",
                                                    varying = list(endogenous),
                                                    timevar = "endogenous",
                                                    v.names = "residual"))

        rownames(residualsL) <- NULL
        residualsL$endogenous <- factor(residualsL$endogenous, levels = 1:n.endogenous, labels = endogenous)
        reorder <- match(interaction(object$sCorrect$old2new.order$XXclusterXX.old, object$sCorrect$old2new.order$XXendogenousXX.old),
                         interaction(residualsL$cluster,residualsL$endogenous))
        return(residualsL[reorder,"residual"])
    }
}

## * residuals.lvmfit2
#' @export
residuals.lvmfit2 <- residuals2.lvmfit2

## * .normalizeResiduals
.normalizeResiduals <- function(epsilon, Omega, type,
                                missing.pattern, unique.pattern, Omega.missing.pattern){
    type <- match.arg(type, choices = c("response","studentized","normalized"), several.ok = FALSE)

    if(type %in% c("studentized")){
        epsilon <- sweep(epsilon,
                         STATS = sqrt(diag(Omega)),
                         FUN = "/",
                         MARGIN = 2)
        ## object$sCorrect$residuals/epsilon
    }else if(type=="normalized"){
        name.endogenous <- colnames(epsilon)
        if(any(is.na(epsilon))==FALSE){
            epsilon <- epsilon %*% matrixPower(Omega, symmetric = TRUE, power = -1/2)
        }else{
            iOmegaHalf.missing.pattern <- lapply(Omega.missing.pattern,matrixPower,symmetric = TRUE, power = -1/2)
            for(iP in names(missing.pattern)){
                iY <- which(unique.pattern[iP,]==1)
                for(iC in missing.pattern[[iP]]){ ## iC <- 1
                    epsilon[iC,iY] <- epsilon[iC,iY] %*% iOmegaHalf.missing.pattern[[iP]]
                }
            }
            
        }
        colnames(epsilon) <- name.endogenous
    }

    return(epsilon)
}

## * .adjustResiduals
.adjustResiduals <- function(epsilon, Psi, Omega, 
                             name.pattern, missing.pattern, unique.pattern,
                             endogenous, n.endogenous, n.cluster){
    if(is.null(Psi)){return(epsilon)}
    n.endogenous <- length(endogenous)

    epsilon.adj <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                          dimnames = list(NULL, endogenous))
    n.pattern <- NROW(unique.pattern)
    
    for(iP in 1:n.pattern){ ## iP <- 1 
        iIndex <- missing.pattern[[iP]]
        iY <- which(unique.pattern[iP,]==1)
        
        iOmega <- Omega[iY,iY,drop=FALSE]
        iPsi <- Psi[iY,iY,drop=FALSE]

        iOmega.chol <- matrixPower(iOmega, symmetric = TRUE, power = 1/2)
        iH <- iOmega %*% iOmega - iOmega.chol %*% iPsi %*% iOmega.chol
        iHM1 <- tryCatch(matrixPower(iH, symmetric = TRUE, power = -1/2), warning = function(w){w})
        if(inherits(iHM1,"warning")){
            stop("Cannot compute the adjusted residuals \n",
                 "Estimated bias too large compared to the estimated variance-covariance matrix \n",
                 "Consider setting argument \'adjust.n\' to FALSE when calling sCorrect \n")
        }
        epsilon.adj[iIndex,iY] <- epsilon[iIndex,iY,drop=FALSE] %*% iOmega.chol %*% iHM1 %*% iOmega.chol
    }

    return(epsilon.adj)
}



######################################################################
### sCorrect-residuals2.R ends here

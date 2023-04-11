### Objective_gaussian_weight.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 17 2020 (16:29) 
## Version: 
## Last-Updated: Jan 12 2022 (12:31) 
##           By: Brice Ozenne
##     Update #: 136
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Estimate LVM With Weights
##' @description Estimate LVM with weights.
##' @name gaussian_weight
##'
##' @param x,object A latent variable model
##' @param data dataset
##' @param estimator name of the estimator to be used
##' @param type must be "cond"
##' @param p parameter value
##' @param weights weight associated to each iid replicate.
##' @param S empirical variance-covariance matrix between variable
##' @param n number of iid replicates
##' @param mu empirical mean
##' @param debug,reindex,mean,constrain,indiv additional arguments not used
##' @param ... passed to lower level functions.
##' 
##' @examples
##' #### linear regression with weights ####
##'
##' ## data
##' df <- data.frame(Y = c(1,2,2,1,2),
##'                  X = c(1,1,2,2,2),
##'                  missing = c(0,0,0,0,1),
##'                  weights = c(1,1,2,1,NA))
##'
##' ## using lm
##' e.lm.GS <- lm(Y~X, data = df)
##' e.lm.test <- lm(Y~X, data = df[df$missing==0,], weights = df[df$missing==0,"weights"])
##' 
##' ## using lvm
##' m <- lvm(Y~X)
##' e.GS <- estimate(m, df)
##' ## e.lava.test <- estimate(m, df[df$missing==0,], weights = df[df$missing==0,"weights"])
##' ## warnings!!
##' e.test <- estimate(m, data = df[df$missing==0,],
##'                    weights = df[df$missing==0,"weights"],
##'                    estimator = "gaussian_weight")
##' 

## * gaussian_weight.estimate.hook
##' @rdname gaussian_weight
##' @export
gaussian_weight.estimate.hook <- function(x, data, estimator, ...){
    dots <- list(...)
    if(identical(estimator,"gaussian_weight")){
        xe <- suppressWarnings(estimate(x, data = data, control = list(iter.max = 0))) ## initialize coefficients
        x$sCorrect <- moments2(xe, data = data, param = NULL, ## initCoef,
                               initialize = TRUE, usefit = FALSE, score = TRUE, information = TRUE, hessian = FALSE, vcov = TRUE, residuals = TRUE, leverage = FALSE, dVcov = FALSE, dVcov.robust = FALSE)
    }
    return( c(list(x=x, data=data, estimator = estimator),dots) )
}

##' @rdname gaussian_weight
##' @export
gaussian_weight_method.lvm <- "nlminb2"

## * gaussian_weight_logLik.lvm
##' @rdname gaussian_weight
##' @export
`gaussian_weight_logLik.lvm` <- function(object, type="cond", p, data, weights,...) {
    ## ** compute mu and Omega
    if(type!="cond"){
        stop("Not implemented for other types than \"cond\"\n ")
    }
    cM <- moments2(object, param = p, data = data, weights = weights,
                   initialize = FALSE, usefit = TRUE, score = FALSE, information = FALSE, hessian = FALSE, vcov = FALSE, residuals = TRUE, leverage = FALSE, dVcov = FALSE, dVcov.robust = FALSE)
    
    ## ** prepare
    name.pattern <- cM$missing$name.pattern
    missing.pattern <- cM$missing$pattern
    unique.pattern <- cM$missing$unique.pattern
    n.pattern <- length(name.pattern)
    
    OmegaM1 <- cM$moment$OmegaM1
    residuals <- cM$residuals
    logLik <- 0

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)
        iResiduals <- residuals[iIndex,iY,drop=FALSE]
        iM <- length(iY)
        if(is.null(weights)){
            logLik <- logLik - (cM$cluster$n.cluster/2) * (iM * log(2*pi) - log(det(iOmegaM1))) - sum((iResiduals %*% iOmegaM1) * iResiduals)/2
        }else{
            logLik <- logLik - (sum(weights)/2) * (iM * log(2*pi) - log(det(iOmegaM1))) - sum(weights[,1]/2 * rowSums((iResiduals %*% iOmegaM1) * iResiduals))
        }
    }
    return(logLik)
}

##' @rdname gaussian_weight
##' @export
`gaussian_weight_objective.lvm` <- function(x, ...) {
    logLik <- gaussian_weight_logLik.lvm(object = x,...)
    return(-logLik)
}

## * gaussian_weight_score.lvm
##' @rdname gaussian_weight
##' @export
gaussian_weight_score.lvm <- function(x, data, p, S, n, mu=NULL, weights=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=TRUE, indiv=FALSE,...) {

    ## if(constrain){
    ##     stop("gaussian_weight_score.lvm does not handle constrain")
    ## }
    ## if(reindex){
    ##     stop("gaussian_weight_score.lvm does not handle reindex")
    ## }
    ## if(!mean){
    ##     stop("gaussian_weight_score.lvm only handles mean")
    ## }
    
    ## ** compute moments
    cM <- moments2(x, param = p, data = data, weights = weights,
                   initialize = FALSE, usefit = TRUE, score = TRUE, information = FALSE, hessian = FALSE, vcov = FALSE, residuals = FALSE, leverage = FALSE, dVcov = FALSE, dVcov.robust = FALSE)

    ## ** export
    if(indiv){
        return(cM$score)
    }else{
        return(colSums(cM$score))
    }
}

## * gaussian_weight_gradient.lvm
##' @rdname gaussian_weight
##' @export
gaussian_weight_gradient.lvm <-  function(...) {
    return(-gaussian_weight_score.lvm(...))
}

## * gaussian_weight_hessian.lvm
##' @rdname gaussian_weight
##' @export
`gaussian_weight_hessian.lvm` <- function(x, p, n, weights=NULL,...) {

    ## ** compute moments
    cM <- moments2(x, param = p, weights = weights,
                   initialize = FALSE, usefit = TRUE, score = FALSE, information = TRUE, hessian = FALSE, vcov = FALSE, residuals = FALSE, leverage = FALSE, dVcov = FALSE, dVcov.robust = FALSE)

    ## ** export
    return(cM$information)
}

######################################################################
### Objective_gaussian_weight.R ends here

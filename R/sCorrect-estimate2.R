### sCorrect-estimate2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: jan 23 2024 (10:30) 
##           By: Brice Ozenne
##     Update #: 2169
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - estimate2
#' @title  Satterthwaite Correction and Small Sample Correction
#' @description Correct the bias of the ML estimate of the variance and compute the first derivative of the information matrix.
#' @name estimate2
#'
#' @param object a \code{lvm} object.
#' @param param [numeric vector, optional] the values of the parameters at which to perform the correction.
#' @param data [data.frame, optional] the dataset relative to which the correction should be performed.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/\code{FALSE}/\code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param tol.max [numeric >0] the largest acceptable absolute difference between two succesive estimates of the bias correction.
#' @param iter.max [integer >0] the maximum number of iterations used to estimate the bias correction.
#' @param derivative [character] should the first derivative of the information matrix be computed using a formula (\code{"analytic"}) or numerical derivative (\code{"numeric"})?
#' @param hessian [logical] should the hessian be stored? Can be \code{NULL} to indicate only if computed during the small sample correction.
#' @param dVcov.robust [logical] should the first derivative of robust variance-covariance matrix be stored?
#' @param trace [logical] should the execution of the function be traced.
#' @param ...  arguments passed to \code{lava::estimate} when using a \code{lvm} object.
#'
#' @details The argument \code{value} is equivalent to the argument \code{bias.correct} of the function \code{summary2}.
#' 
#' @concept estimator
#' @keywords smallSampleCorrection
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' 
#' #### latent variable model ####
#' m.lvm <- lvm(Y1~X1+X2+Z1)
#'
#' e2.lvm <- estimate2(m.lvm, data = dW)
#' summary2(e2.lvm)
#' 
#' @export
`estimate2` <-
    function(object, param, data,
             ssc, df,
             derivative, hessian, dVcov.robust,
             iter.max, tol.max, trace , ...) UseMethod("estimate2")


## * estimate2.lvm
#' @rdname estimate2
#' @export
estimate2.lvm <- function(object, param = NULL, data = NULL,
                          ssc = lava.options()$ssc, df = lava.options()$df,
                          derivative = "analytic", hessian = FALSE, dVcov.robust = FALSE,
                          iter.max = 100, tol.max = 1e-6, trace = 0, ...){

    if(!is.null(param)){
        warning("Argument \'param\' not used with lvm objects. \n")
    }
    
    out <- lava::estimate(x = object, data = data, ...)
    return(estimate2(out, param = NULL, data = NULL,
                     ssc = ssc, df = df,
                     derivative = derivative, hessian = hessian, dVcov.robust = FALSE,
                     iter.max = iter.max, tol.max = tol.max, trace = trace))

}

## * estimate2.lvmfit
#' @rdname estimate2
#' @export
estimate2.lvmfit <- function(object, param = NULL, data = NULL,
                             ssc = lava.options()$ssc, df = lava.options()$df,
                             derivative = "analytic", hessian = FALSE, dVcov.robust = FALSE,
                             iter.max = 100, tol.max = 1e-6, trace = 0, ...){

    ## ** preliminary tests
    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    if("multigroupfit" %in% class(object)){
        stop("estimate2 cannot handle multigroup models \n")
    }

    if(inherits(object,"lvmfit") && length(object$model$attributes$ordinal)>0){
        name.t <- names(object$model$attributes$ordinal)
        stop("estimate2 does not handle ordinal variables \n",
             "ordinal variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    if(inherits(object,"lvmfit") && length(object$model$attributes$transform)>0){
        name.t <- names(object$model$attributes$transform)
        stop("estimate2 does not handle transformed variables \n",
             "transformed variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }

    ## arguments
    if(identical(ssc,FALSE) || identical(ssc,NA)){
        ssc <- "none"
    }
    if(identical(df,FALSE) || identical(df,NA)){
        df <- "none"
    }
    ssc <- match.arg(tolower(ssc), c("none","residuals","residuals0","cox"))
    df <- match.arg(tolower(df), c("none","satterthwaite"))

    if(df %in% "satterthwaite" || ssc %in% "cox"){
        second.order <- TRUE
    }else if(ssc %in% c("residuals","residuals0")){
        second.order <- FALSE
    }else{
        second.order <- FALSE
    }

    ## ** initialize object
    if(trace>0){cat("Initialization:")}
    object$sCorrect <- moments2(object, data = data, param = param, Psi = NULL,
                                initialize = TRUE, usefit = TRUE,
                                score = TRUE, information = TRUE, hessian = hessian, vcov = TRUE,
                                dVcov = (ssc == "cox")  || (ssc == "none" && df == "satterthwaite"), dVcov.robust = (ssc == "none" && df == "satterthwaite" && dVcov.robust),
                                residuals = TRUE, leverage = FALSE, derivative = derivative)  ## setting leverage to FALSE is like initialization to 0
    if(trace>0){cat(" done \n")}

    ## ** bias correction    
    if(ssc != "none"){

        ## *** initialize bias correction
        if(trace>0){cat("Initialize bias correction \n")}
        if(ssc=="Cox"){
            object.ssc <- list(type = "Cox",
                               param0 = object$sCorrect$param,
                               Omega0 = object$sCorrect$moment$Omega)
        }else if(ssc %in% c("residuals","residuals0")){
            object.ssc <- .init_sscResiduals(object)
        }
        
        ## *** perform bias correction
        if(trace>0){cat("Perform bias correction \n")}
        iCV <- FALSE
        iIter <- 0
        iTol <- Inf
        iiParam <- object$sCorrect$param
        ## cat(iTol," (",iiParam,") \n")
            
        while(iCV == FALSE && iIter < iter.max){
            if(trace>0){cat("*")}

            ## bias correction
            if(ssc == "Cox"){
                iParam <- .sscCoxSnell(object, ssc = object.ssc)
                object.ssc$JJK <- attr(iParam,"JJK")
                object.ssc$lm <- attr(iParam,"lm")
                object.ssc$Psi <- object$sCorrect$moment$Omega - object.ssc$Omega0
            }else if(ssc %in% c("residuals","residuals0")){
                iParam <- .sscResiduals(object, ssc = object.ssc)
                object.ssc$Omega <- attr(iParam,"Omega")
                object.ssc$Psi <- attr(iParam,"Psi")
                ## use previous Omega to compute leverage and residuals
                attr(object.ssc$Omega,"Omega.leverage") <- object$sCorrect$moment$Omega
                attr(object.ssc$Omega,"dOmega.leverage") <- object$sCorrect$dmoment$dOmega
                attr(object.ssc$Omega,"Omega.residuals") <- object$sCorrect$moment$Omega
            }
            ## object.ssc$Omega0 + object.ssc$Psi - object.ssc$Omega

            ## cv criteria
            iIter <- iIter + 1
            iTol <- max(abs(iParam-iiParam))
            ## cat(iTol," (",iParam,") \n")
            iiParam <- iParam
            iCV <- iTol <= tol.max

            ## update moments
            ## if ssc=="residuals0" then do not rescale the residuals according the the bias
            if(iCV==FALSE && iIter < iter.max){
                object$sCorrect <- moments2(object, param = iParam, Psi = if(ssc!="residuals0"){object.ssc$Psi}else{NULL}, Omega = object.ssc$Omega, 
                                            initialize = FALSE, usefit = TRUE,
                                            score = TRUE, information = TRUE, hessian = FALSE, vcov = TRUE,
                                            dVcov = (ssc == "cox"), dVcov.robust = FALSE,
                                            residuals = TRUE, leverage = TRUE, derivative = derivative)
            }else{
                object$sCorrect <- moments2(object, param = iParam, Psi = if(ssc!="residuals0"){object.ssc$Psi}else{NULL}, Omega = object.ssc$Omega, 
                                            initialize = FALSE, usefit = TRUE,
                                            score = TRUE, information = TRUE, hessian = hessian, vcov = TRUE,
                                            dVcov = (df == "satterthwaite"), dVcov.robust = dVcov.robust,
                                            residuals = TRUE, leverage = TRUE, derivative = derivative)

                object$sCorrect$ssc <- c(object.ssc,
                                         cv = iCV,
                                         iter = iIter,
                                         tol = iTol,
                                         iter.max = iter.max,
                                         tol.max = tol.max
                                         )
                
            }
        }
        
        ## *** assess convergence
        if(iCV == FALSE){
            warning("small sample correction did not reach convergence after ",iIter," iteration",if(iIter>1){"s"}else{""},". \n")
        }
        if(trace > 0){
            cat("\n")
        }

    }else{
        object$sCorrect$ssc$type <- "none"
    }

    ## ** degrees of freedom    
    object$sCorrect$df <- df ## degrees of freedom are computed later (by compare2)
    
    ## ** restaure original param order
    name.param <- object$sCorrect$name.param
    if(!is.null(name.param)){
        object$sCorrect$param <- object$sCorrect$param[name.param]
        names(object$sCorrect$param) <- names(name.param)

        if(!is.null(object$sCorrect$score)){
            object$sCorrect$score <- object$sCorrect$score[,name.param,drop=FALSE]
            colnames(object$sCorrect$score) <- names(name.param)
        }
        if(!is.null(object$sCorrect$information)){
            object$sCorrect$information <- object$sCorrect$information[name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$information) <- list(names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$vcov.param)){
            object$sCorrect$vcov.param <- object$sCorrect$vcov.param[name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$vcov.param) <- list(names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$hessian)){
            object$sCorrect$hessian <- object$sCorrect$hessian[name.param,name.param,,drop=FALSE]
            dimnames(object$sCorrect$hessian) <- list(names(name.param),names(name.param),NULL)
        }
        if(!is.null(object$sCorrect$dInformation)){
            object$sCorrect$dInformation <- object$sCorrect$dInformation[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dInformation) <- list(names(name.param),names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$dVcov.param)){
            object$sCorrect$dVcov.param <- object$sCorrect$dVcov.param[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dVcov.param) <- list(names(name.param),names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$dRvcov.param)){
            object$sCorrect$dRvcov.param <- object$sCorrect$dRvcov.param[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dRvcov.param) <- list(names(name.param),names(name.param),names(name.param))
        }
    }
    
    ## ** export
    class(object) <- append("lvmfit2",class(object))
    return(object)    
}

## * estimate2.list
#' @rdname estimate2
#' @export
estimate2.list <- function(object, ...){
    object.class <- class(object)
    object <- lapply(object, estimate2, ...)
    class(object) <- object.class
    return(object)
}

## * estimate2.mmm
#' @rdname estimate2
#' @export
estimate2.mmm <- estimate2.list

##----------------------------------------------------------------------
### sCorrect.R ends here










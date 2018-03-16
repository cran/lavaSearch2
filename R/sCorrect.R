### sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: mar 15 2018 (18:15) 
##           By: Brice Ozenne
##     Update #: 931
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - sCorrect
#' @title  Compute the Derivative of the Information Matrix
#' @description Compute the derivative of the information matrix.
#' @name sCorrect
#'
#' @param object,x a \code{gls}, \code{lme}, or \code{lvm} object.
#' @param param [numeric vector, optional] the values of the parameters at which to perform the correction.
#' @param data [data.frame, optional] the dataset relative to which the correction should be performed.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param value [logical] value for the arguments \code{adjust.Omega} and \code{adjust.n}.
#' @param df [logical] should the first derivative of the expected information matrix be computed. Required when computing the degrees of freedom of the test statistics.
#' @param adjust.Omega [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param adjust.n [logical] should the correction for the degree of freedom be performed?
#' @param tol [numeric >0] the minimum absolute difference between two estimation of the small sample bias.
#' Below this value, the algorithm used to estimate the bias stop.
#' @param n.iter [integer >0] the maximum number of iterations used to estimate the small sample bias of the residual variance-covariance matrix. 
#' @param numeric.derivative [logical] should a numerical derivative be used to compute the first derivative of the information matrix?
#' Otherwise an analytic formula is used.
#' @param trace [logical] should the execution of the function be traced.
#' @param score [internal] export the score.
#' @param ... [internal] only used by the generic method or by the <- methods.
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#' 
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' system.time(
#' sCorrect(e.lm) <- TRUE
#')
#' 
#' ## gls model
#' library(nlme)
#' e.gls <- gls(formula.lvm, data = d, method = "ML")
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE
#' summary2(e.gls)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' sCorrect(e.lvm) <- TRUE
#' summary2(e.lvm)
#' 
#' @export
`sCorrect` <-
  function(object, ...) UseMethod("sCorrect")


## * sCorrect.lm
#' @rdname sCorrect
#' @export
sCorrect.lm <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                        score = TRUE, df = TRUE, numeric.derivative = FALSE,
                        param = NULL, data = NULL,
                        tol = 1e-5, n.iter = 20, trace = 0, ...){
    
### ** Extract quantities from object
    name.endogenous <- all.vars(stats::update(formula(object), ".~1"))
    n.cluster <- stats::nobs(object) + length(object$na.action)

    if(is.null(param)){
        param <- .coef2(object)
        param["sigma2"] <- mean(residuals(object)^2)
        model.param <- param
    }else{
        model.param <- .coef2(object)
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }
    name.param <- names(model.param)
    name.meanparam <- attr(model.param,"mean.coef")
    name.varparam <- attr(model.param,"var.coef")

### ** Compute conditional moments
    if(trace>0){
        cat("Compute conditional moments")
    }
    dMoments <- conditionalMoment(object, name.endogenous = name.endogenous,
                                  second.order = df)
    if(trace>0){
        cat(" - done \n")
    }
    
    ### ** Compute residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    if(is.null(data)){
        X <- model.matrix(object)
    }else{
        X <- model.matrix(formula(object), data)
    }
    Y <- object$residuals + object$fitted.values
    object.residuals <- Y - X %*% cbind(model.param[attr(model.param,"mean.coef")])    
    dimnames(object.residuals) <- list(NULL, name.endogenous)
    if(trace>0){
        cat("- done \n")
    }
    
    ## ** Compute residual variance
    if(trace>0){
        cat("* Reconstruct estimated residual variance ")
    }
    Omega <- matrix(param["sigma2"], nrow = 1, ncol = 1,
                    dimnames = list(name.endogenous, name.endogenous))
    if(trace>0){
        cat("- done \n")
    }

    ## ** args
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,
                 df = df,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter)
    
    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

    out <- .sCorrect(object,
                     data = data,
                     param = model.param,
                     epsilon = object.residuals,
                     Omega = Omega,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = NULL,
                     d2Omega = NULL,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = dMoments$name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = NULL,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     args = args,
                     trace = trace,
                     ...)
    
    ## ** export
    return(out)    
}

## * sCorrect.gls
#' @rdname sCorrect
#' @export
sCorrect.gls <- function(object, cluster, adjust.Omega = TRUE, adjust.n = TRUE,
                         score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                         param = NULL, data = NULL,
                         tol = 1e-5, n.iter = 20, trace = 0,
                         ...){

### ** limitations
    if(object$method!="ML"){
        if(adjust.Omega==TRUE || adjust.n == TRUE){
            warning("Small sample corrections were derived for ML not for REML\n")
        }else{
            warning("The Satterthwaite approximation ignores that fact that the model was fitted using REML\n")
        }
    }
    
    ## check valid class for corStruct and varStruct: see .getVarCov2
### ** Extract quantities from the model
    ## *** data
    if(is.null(data)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()))
    }
    
    ## *** endogenous variable
    formula.object <- .getFormula2(object)
    name.Y <- all.vars(stats::update(formula.object, ".~1"))
    Y <- data[[name.Y]]
    
    ## *** parameters
    model.param <- .coef2(object)
    if(!is.null(param)){        
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }
    name.param <- names(model.param)
    name.meanparam <- attr(model.param,"mean.coef")
    name.varparam <- c(attr(model.param,"var.coef"),
                       attr(model.param,"cor.coef"),
                       attr(model.param,"ran.coef"))

    ## *** group
    if(trace>0){
        cat("* Reconstruct iid sample ")
    }
    res.cluster <- .getCluster2(object,
                                data = data,
                                cluster = cluster)
    if(trace>0){
        cat("- done \n")
    }
    
    ## *** repetition relative to each observation
    if(trace>0){
        cat("* Relate observations to endogenous variables ")
    }

    res.index <- .getIndexOmega2(object,
                                 param = model.param,
                                 attr.param = attributes(model.param),
                                 name.Y = name.Y,
                                 cluster = res.cluster$cluster,
                                 levels.cluster = res.cluster$levels.cluster,
                                 data = data)
    index.Omega <- res.index$index.Omega

    if(trace>0){
        cat("- done \n")
    }

    ### ** Reconstruct residuals variance covariance matrix
    if(trace>0){
        cat("* Reconstruct estimated residual variance-covariance matrix ")
    }
    
    Omega <- .getVarCov2(object,
                         param = model.param,
                         attr.param = attributes(model.param),
                         name.endogenous = res.index$name.endogenous,
                         n.endogenous = res.index$n.endogenous,
                         ref.group = res.index$ref.group)

    if(trace>0){
        cat("- done \n")
    }
    cluster <- res.cluster$cluster
    n.cluster <- res.cluster$n.cluster
    name.endogenous <- res.index$name.endogenous
    n.endogenous <- res.index$n.endogenous

### ** Compute conditional moments and derivatives
    if(trace>0){
        cat("* Compute conditional moments ")
    }
    dMoments <- conditionalMoment(object,
                                  data = data,
                                  formula = formula.object,
                                  param = model.param,
                                  attr.param = attributes(model.param)[-1],
                                  ref.group = res.index$ref.group,
                                  second.order = df,
                                  n.cluster = n.cluster,
                                  cluster = cluster,
                                  name.endogenous = name.endogenous,
                                  index.Omega = index.Omega)
    if(trace>0){
        cat("- done \n")
    }
    
### ** Compute observed residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))
    for(iC in 1:n.cluster){
        iIndex <- which(cluster == iC)
        epsilon[iC,index.Omega[[iC]]] <- as.double(Y[iIndex] - dMoments$X[iIndex,,drop=FALSE] %*% model.param[name.meanparam])
    }

    if(trace>0){
        cat("- done \n")
    }
    ## stats::residuals(object)-as.double(t(epsilon))

    ## ** Check missing value
    if(all(!is.na(epsilon))){
        index.Omega <- NULL
    }
    
    ## ** param with non-zero third derivative
    name.3deriv <- name.varparam

    ## ** args
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,                     
                 df = df,
                 cluster = cluster,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter)
    
    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }
    out <- .sCorrect(object,
                     data = data,
                     param = model.param,
                     epsilon = epsilon,
                     Omega = Omega,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = NULL,
                     d2Omega = dMoments$d2Omega,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     args = args,
                     trace = trace,
                     ...)
    
    ## ** export    
    return(out)          
 
}

## * sCorrect.lme
#' @rdname sCorrect
#' @export
sCorrect.lme <- sCorrect.gls

## * sCorrect.lvmfit
#' @rdname sCorrect
#' @export
sCorrect.lvmfit <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                            score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                            param = NULL, data = NULL,
                            tol = 1e-5, n.iter = 20, trace = 0,
                            ...){

    ## ** Check valid lvm object
    if("multigroupfit" %in% class(object)){
        stop("sCorrect does not handle multigroup models \n")
    }
   
    if(length(object$model$attributes$ordinal)>0){
        name.t <- names(object$model$attributes$ordinal)
        stop("sCorrect does not handle ordinal variables \n",
             "ordinal variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    if(length(object$model$attributes$transform)>0){
        name.t <- names(object$model$attributes$transform)
        stop("sCorrect does not handle transformed variables \n",
             "transformed variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    ## ** Extract quantities from object
    name.endogenous <- endogenous(object)

    model.param <- lava::pars(object)
    if(!is.null(param)){
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }

    if(is.null(data)){
        data <- as.data.frame(object$data$model.frame)
    }

    name.param <- names(model.param)

    n.latent <- length(latent(object))

    ### ** number of samples
    test.NNA <- rowSums(is.na(data[,name.endogenous,drop=FALSE]))==0    
    if(any(test.NNA==FALSE) && !inherits(object,"lvm.missing")){ ## complete case analysis
        if(trace>0){
            cat("* Exclude missing values and recompute moments and residuals ")
        }        
        data <- data[which(test.NNA),,drop=FALSE]
        if(trace>0){
            cat("- done \n")
        }        
    }
    
    n.cluster <- NROW(data)

    ### ** Compute conditional moments and derivatives
    if(trace>0){
        cat("* Compute conditional moments ")
    }
    dMoments <- conditionalMoment(object, data = data, param = model.param,
                                  second.order = df, usefit = TRUE)
    if(trace>0){
        cat("- done \n")
    }

    name.meanparam <- names(dMoments$dmu)
    name.varparam <- names(dMoments$dOmega)

    #### ** Compute residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    epsilon <- .calcResidualsLVM(data = data, dMoments = dMoments,
                                 n.latent = n.latent,
                                 name.endogenous = name.endogenous)
    if(trace>0){
        cat("- done \n")
    }

### ** Identify missing values
    if(any(test.NNA==FALSE) && inherits(object,"lvm.missing")){ ## full information
            if(trace>0){
                cat("* Identify missing values ")
            }
            index.Omega <- lapply(1:n.cluster,function(iC){which(!is.na(epsilon[iC,]))})
            if(trace>0){
                cat("- done \n")
            }        
    }else{
        index.Omega <- NULL
    }

    ### ** Compute residual variance covariance matrix
    if(trace>0){
        cat("* Reconstruct residual variance-covariance matrix ")
    }

    Omega <- .calcOmegaLVM(dMoments, n.latent = n.latent)

    if(trace>0){
        cat("- done \n")
    }

    ## ** param with non-zero third derivative
    type.3deriv <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
    index.keep <- intersect(which(!is.na(dMoments$df.param$lava)),
                            which(dMoments$df.param$detail %in% type.3deriv)
                            )
    
    name.3deriv <- dMoments$df.param[index.keep, "originalLink"]

    ## ** args
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,                     
                 df = df,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter)

    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

   out <- .sCorrect(object,
                     data = data,
                     param = model.param,
                     epsilon = epsilon,
                     Omega = Omega,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = dMoments$d2mu,
                     d2Omega = dMoments$d2Omega,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     args = args,
                     trace = trace,
                     ...)

    ## ** export
    return(out)       
}

## * sCorrect.lvmfit2
#' @rdname sCorrect
#' @export
sCorrect.lvmfit2 <- function(object, ...){
    class(object) <- setdiff(class(object),"lvmfit2")
    return(sCorrect(object, ...))    
}
## * .sCorrect
.sCorrect <- function(object, data, param, epsilon, Omega, dmu, dOmega, d2mu, d2Omega, 
                      name.param, name.meanparam, name.varparam, name.endogenous, name.3deriv,
                      n.cluster, index.Omega,
                      adjust.Omega, adjust.n, tol, n.iter, score, derivative, args, trace){

    n.param <- length(param)
    if(!is.null(index.Omega)){
        n.endogenous.cluster <- lapply(index.Omega,length)        
    }else{
        n.endogenous.cluster <- NULL
    }
    ## ** order param names
    name.meanparam <- as.character(sort(factor(name.meanparam, levels = name.param)))
    if(any(is.na(name.meanparam))){
        stop("An element in name.meanparam is not in name.param. \n")
    }
    if(length(name.meanparam)>0 && !identical(sort(name.meanparam),sort(names(dmu)))){
        stop("Mismatch first derivative of the conditional mean and name.meanparam \n")
    }
    dmu <- dmu[name.meanparam]    

    name.varparam <- as.character(sort(factor(name.varparam, levels = name.param)))
    if(any(is.na(name.varparam))){
        stop("An element in name.varparam is not in name.param. \n")
    }
    if(length(name.varparam)>0 && !identical(sort(name.varparam),sort(names(dOmega)))){
        stop("Mismatch first derivative of the conditional variance and name.varparam \n")
    }
    dOmega <- dOmega[name.varparam]

    ## ** corrected ML estimates
    out  <- adjustEstimate(epsilon = epsilon,
                           Omega = Omega,
                           dmu = dmu,
                           dOmega = dOmega,
                           n.cluster = n.cluster,
                           name.param = name.param,
                           name.meanparam = name.meanparam,
                           name.varparam = name.varparam,
                           name.endogenous = name.endogenous,
                           index.Omega = index.Omega, ## mode2
                           adjust.Omega = adjust.Omega,
                           adjust.n = adjust.n,
                           tol = tol, n.iter = n.iter,
                           trace = trace)    
    out$param <- param

    ## ** corrected score
    if(score){
        if(trace>0){
            if(adjust.n == FALSE && adjust.Omega == FALSE){
                cat("* Compute score ")
            }else{
                cat("* Compute corrected score ")
            }
        }
        out$score <- .score2(epsilon = out$epsilon,
                             Omega = out$Omega,
                             OmegaM1 = out$OmegaM1,
                             dmu = dmu,
                             dOmega = dOmega,
                             name.param = name.param,
                             name.meanparam = name.meanparam,
                             name.varparam = name.varparam,
                             index.Omega = index.Omega, ## mode2
                             n.cluster = n.cluster,
                             indiv = TRUE)
        if(trace>0){
            cat("- done \n")
        }
    }
    
    ## ** first derivative of the expected information matrix
    if(length(name.3deriv)==0){
        out$dVcov.param <- NULL
    }else if(derivative == "none"){
        out$dVcov.param <- NA
    }else if(derivative == "numeric"){
        if(trace>0){
            cat("Compute first derivative of the information matrix using numerical differentiation ")
        }
        if(adjust.Omega || adjust.n){
            warning("The numerical derivative of the information matrix is computed ignoring the small sample correction \n")
        }

        if("lvmfit" %in% class(object)){
            object$conditionalMoment <- conditionalMoment(lava::Model(object), data = data,
                                                          usefit = FALSE, second.order = FALSE,
                                                          name.endogenous = endogenous(object),
                                                          name.latent = latent(object))
        }
        args.tempo <- args
        args.tempo$data <- data
        args.tempo$df <- FALSE
        args.tempo$score <- FALSE

        ## *** direct computation of the variance-covariance matrix
        calcVcov <- function(iParam){ # x <- p.obj
            pp <- param
            pp[names(iParam)] <- iParam
            
            vcov.param <- do.call(sCorrect,
                                  args = c(list(object, param = pp), args.tempo))$vcov.param
            return(vcov.param)
        }

        ### *** numerical derivative
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }
        jac.param <- param[name.3deriv]
        res.numDeriv <- numDeriv::jacobian(calcVcov, x = jac.param, method = "Richardson")
        
        out$dVcov.param <- array(res.numDeriv,
                                 dim = c(n.param,n.param,length(name.3deriv)),
                                 dimnames = list(name.param, name.param, name.3deriv))
         if(trace>0){
            cat("- done \n")
         }
        
    }else if(derivative == "analytic"){
        if(trace>0){
            cat("* Compute first derivative of the information matrix using analytic formula ")
        }

        dInfo.dtheta <- .dInformation2(dmu = dmu,
                                       d2mu = d2mu,
                                       dOmega = dOmega,
                                       d2Omega = d2Omega,
                                       Omega = out$Omega,
                                       OmegaM1 = out$OmegaM1,
                                       n.corrected = out$n.corrected,
                                       n.cluster = n.cluster,
                                       index.Omega = index.Omega,
                                       leverage = out$leverage,
                                       name.param  = name.param,
                                       name.3deriv = name.3deriv)
        p3 <- dim(dInfo.dtheta)[3]
        out$dVcov.param <- array(NA, dim = dim(dInfo.dtheta), dimnames = dimnames(dInfo.dtheta))
        
        for(iP in 1:p3){
            out$dVcov.param[,,iP] <- - out$vcov.param %*% dInfo.dtheta[,,iP] %*% out$vcov.param 
        }

        if(trace>0){
            cat("- done \n")
        }
    }
       
    ## ** export
    out$args <- args
    return(out)
}

## * sCorrect<-
#' @rdname sCorrect
#' @export
`sCorrect<-` <-
  function(x, ..., value) UseMethod("sCorrect<-")

## * sCorrect<-.lm
#' @rdname sCorrect
#' @export
`sCorrect<-.lm` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lm2",class(x))

    return(x)
}    
## * sCorrect<-.gls
#' @rdname sCorrect
#' @export
`sCorrect<-.gls` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("gls2",class(x))

    return(x)
}    
## * sCorrect<-.lme
#' @rdname sCorrect
#' @export
`sCorrect<-.lme` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lme2",class(x))
    
    return(x)
}    

## * sCorrect<-.lvmfit
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit` <- function(x, ..., value){
    x$sCorrect <-  try(sCorrect(x, ..., adjust.Omega = value, adjust.n = value), silent = TRUE)
    if(value == TRUE && inherits(x$sCorrect,"try-error")){
        warn <- x$sCorrect
        x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = FALSE)
        attr(x$sCorrect,"warning") <- warn
        warning("sCorrect failed and has been re-run setting the argument \'adjust.n\' to FALSE \n",
                "see the attribute \"warning\" of object$sCorrect for the error message \n")
    }
    class(x) <- append("lvmfit2",class(x))

    return(x)
}    

## * sCorrect<-.lvmfit2
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit2` <- function(x, ..., value){

    class(x) <- setdiff(class(x),"lvmfit2")
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lvmfit2",class(x))
    
    return(x)
}

## * .calcResidualsLVM
.calcResidualsLVM <- function(data, dMoments, n.latent, name.endogenous){

    ## ** fitted value
    object.fitted <- dMoments$value$nu.XK
    if(n.latent>0){
        object.fitted <- object.fitted + dMoments$value$alpha.XGamma.iIB %*% dMoments$value$Lambda
    }

    ## ** residuals
    out <- data[, name.endogenous] - object.fitted
    
    ## ** export
    return(as.matrix(out))

}

## * .calcOmegaLVM
.calcOmegaLVM <- function(dMoments, n.latent){
    if(n.latent>0){
        Omega <- dMoments$value$tLambda.tiIB.Psi.iIB %*% dMoments$value$Lambda + dMoments$value$Sigma
    }else{
        Omega <- dMoments$value$Sigma
    }
    return(Omega)
}

##----------------------------------------------------------------------
### sCorrect.R ends here








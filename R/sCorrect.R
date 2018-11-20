### sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: nov 11 2018 (14:32) 
##           By: Brice Ozenne
##     Update #: 1350
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - sCorrect
#' @title  Satterthwaite Correction and Small Sample Correction
#' @description Correct the bias of the ML estimate of the variance and compute the first derivative of the information matrix.
#' @name sCorrect
#'
#' @param object,x a \code{gls}, \code{lme}, or \code{lvm} object.
#' @param param [numeric vector, optional] the values of the parameters at which to perform the correction.
#' @param data [data.frame, optional] the dataset relative to which the correction should be performed.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param value [logical] value for the arguments \code{adjust.Omega} and \code{adjust.n}.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
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
#' @details The argument \code{value} is equivalent to the argument \code{bias.correct} of the function \code{summary2}.
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
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' system.time(
#' sCorrect(e.lm) <- TRUE ## i.e. bias.correct = TRUE
#')
#' 
#' ## gls model
#' library(nlme)
#' e.gls <- gls(formula.lvm, data = d, method = "ML")
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE ## i.e. bias.correct = TRUE
#' summary2(e.gls)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' sCorrect(e.lvm) <- TRUE ## i.e. bias.correct = TRUE
#' summary2(e.lvm)
#' 
#' @export
`sCorrect` <-
    function(object, adjust.Omega, adjust.n,
             score, df, numeric.derivative,
             param, data,
             tol, n.iter, trace, ...) UseMethod("sCorrect")


## * sCorrect.lm
#' @rdname sCorrect
#' @export
sCorrect.lm <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                        score = TRUE, df = TRUE, numeric.derivative = FALSE,
                        param = NULL, data = NULL,
                        tol = 1e-5, n.iter = 20, trace = 0, ...){
    
### ** Extract quantities from object
    name.endogenous <- all.vars(stats::update(formula(object), ".~1"))

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

    if(is.null(data)){
      data <- model.frame(object)
    }

### ** Number of samples
    test.NNA <- sum(is.na(data[[name.endogenous]]))==0
    if(any(test.NNA==FALSE)){ ## complete case analysis
        if(trace>0){
            cat("* Exclude missing values and recompute moments and residuals ")
        }        
        data <- data[which(test.NNA),,drop=FALSE]
        if(trace>0){
            cat("- done \n")
        }        
    }
    
    n.cluster <- NROW(data)
    
### ** Compute conditional moments
    if(trace>0){
        cat("Compute conditional moments")
    }
    object$conditionalMoment <- conditionalMoment(object, data = data, param = model.param,
                                                  name.endogenous = name.endogenous,
                                                  first.order = TRUE, second.order = FALSE)
    if(trace>0){
        cat(" - done \n")
    }
    
    ### ** Compute residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    object.residuals <- data[[name.endogenous]] - object$conditionalMoment$mu    
    dimnames(object.residuals) <- list(NULL, name.endogenous)
    if(trace>0){
        cat("- done \n")
    }

### ** args
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,
                 df = df,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter)

    if(df && numeric.derivative){
        argsNumDeriv <- list(data=data)
    }else{
        argsNumDeriv <- list()
    }
    
### ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = object.residuals,
                     name.param = name.param,
                     name.endogenous = name.endogenous,
                     n.cluster = n.cluster,
                     index.Omega = NULL,
                     score = score,
                     derivative = derivative,
                     args = args,
                     argsNumDeriv = argsNumDeriv,
                     trace = trace,
                     ...)
    
    ## ** export
    return(out)    
}

## * sCorrect.lm2
#' @rdname sCorrect
#' @export
sCorrect.lm2 <- function(object, ...){
    class(object) <- setdiff(class(object),"lm2")
    return(sCorrect(object, ...))    
}
## * sCorrect.gls
#' @rdname sCorrect
#' @export
sCorrect.gls <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                         score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                         param = NULL, data = NULL,
                         tol = 1e-5, n.iter = 20, trace = 0,
                         cluster, ...){
### ** limitations
    if(object$method!="ML"){
        if(adjust.Omega==TRUE || adjust.n == TRUE){
            warning("Small sample corrections were derived for ML not for REML\n")
        }else if(df){
            warning("The Satterthwaite approximation ignores that fact that the model was fitted using REML\n")
        }
    }
    
    ## check valid class for corStruct and varStruct: see .getVarCov2
### ** Extract quantities from the model

    ## *** data
    if(is.null(data)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()))
        
        if(length(object$na.action)>0){ ## remove rows corresponding to missing values
            data <- data[setdiff(1:NROW(data),object$na.action),,drop=FALSE]
        }
    }
    
    ## *** endogenous variable
    formula.object <- .getFormula2(object)
    name.Y <- all.vars(stats::update(formula.object, ".~1"))
    
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
    n.cluster <- res.cluster$n.cluster
    cluster <- res.cluster$cluster
    if(length(cluster) != NROW(data)){
        if(length(object$na.action)>0){
            stop("Number of rows of \'data\' does not match length of cluster \n",
                 "Consider removing the rows of \'data\' containing NA before fitting the model \n")
        }else{
            stop("Number of rows of data does not match length of cluster \n")
        }
    }
    if(trace>0){
        cat("- done \n")
    }

    ## data before re-ordering
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,                     
                 df = df,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter,
                 cluster = cluster) ## for score2

    if(df && numeric.derivative){
        argsNumDeriv <- list(data = data)
    }else{
        argsNumDeriv <- list()
    }
    
    ## *** repetition relative to each observation
    if(trace>0){
        cat("* Relate observations to endogenous variables ")
    }
    res.index <- .getIndexOmega2(object,
                                 param = model.param,
                                 attr.param = attributes(model.param),
                                 name.Y = name.Y,
                                 cluster = cluster,
                                 levels.cluster = res.cluster$levels.cluster,
                                 data = data)
    name.endogenous <- res.index$name.endogenous
    n.endogenous <- res.index$n.endogenous
    index.Omega <- res.index$index.Omega
    if(trace>0){
        cat("- done \n")
    }
    
    ## *** sort data by group
    vec.endogenous <- rep(NA, length(cluster))
    for(iC in 1:n.cluster){ ## iC <- 1        
        ## vec.endogenous[cluster==res.cluster$levels.cluster[iC]] <- index.Omega[[res.cluster$levels.cluster[iC]]]
        vec.endogenous[cluster==iC] <- index.Omega[[iC]]
    }
    order.obs <- order(cluster,vec.endogenous)
    if(is.unsorted(order.obs)==TRUE){
        test.reorder <- TRUE
        data <- data[order.obs,,drop=FALSE]

        cluster <- cluster[order.obs]
        vec.endogenous <- vec.endogenous[order.obs]
        res.cluster$levels.cluster <- unique(cluster)        
    }else{
        test.reorder <- FALSE
    }
    ## for vector format to matrix format (for residuals and fitted values)
    vec.OmegaMat <- cluster + (vec.endogenous-1)*n.cluster
    ## M.check <- matrix(NA, nrow = n.cluster, ncol = n.endogenous)
    ## M.check[vec.endogenous + (cluster-1)*n.endogenous] <- data[["G"]]
    ## M.check[cluster + (vec.endogenous-1)*n.cluster] <- data[["G"]]

### ** Compute conditional moments and derivatives
    if(trace>0){
        cat("* Compute conditional moments ")
    }
    object$conditionalMoment <- conditionalMoment(object,
                                                  data = data,
                                                  formula = formula.object,
                                                  param = model.param,
                                                  attr.param = attributes(model.param)[-1],
                                                  ref.group = res.index$ref.group,
                                                  first.order = TRUE,
                                                  second.order = FALSE,
                                                  n.cluster = n.cluster,
                                                  cluster = cluster,
                                                  name.endogenous = name.endogenous,
                                                  n.endogenous = n.endogenous,
                                                  index.Omega = index.Omega,
                                                  vec.OmegaMat = vec.OmegaMat)
    if(trace>0){
        cat("- done \n")
    }
   
### ** Compute observed residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))
    epsilon[vec.OmegaMat] <- data[[name.Y]]
    epsilon <- epsilon -  object$conditionalMoment$mu
    
    if(trace>0){
        cat("- done \n")
    }
    ## stats::residuals(object)-na.omit(as.double(t(epsilon)))
    ## epsilon + object$conditionalMoment$mu
    ## data
    ## as.double(stats::residuals(object))
    ## ** Check missing value
    if(all(!is.na(epsilon))){
        index.Omega <- NULL
    }
    
    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }
    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = epsilon,
                     name.param = name.param,
                     name.endogenous = name.endogenous,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     score = score,
                     derivative = derivative,
                     args = args,
                     argsNumDeriv = argsNumDeriv,
                     trace = trace,
                     ...)
    
    ## ** export
    ## *** restaure original order
    if(test.reorder==TRUE){
        restaure.order <- order(order.obs[!duplicated(cluster)])
        out$score <- out$score[restaure.order,,drop=FALSE]
        out$residuals <- out$residuals[restaure.order,,drop=FALSE]
        out$leverage <- out$leverage[restaure.order,,drop=FALSE]
    }    
    ##
    return(out)          
 
}

## * sCorrect.gls2
#' @rdname sCorrect
#' @export
sCorrect.gls2 <- function(object, ...){
    class(object) <- setdiff(class(object),"gls2")
    return(sCorrect(object, ...))    
}
## * sCorrect.lme
#' @rdname sCorrect
#' @export
sCorrect.lme <- sCorrect.gls
## * sCorrect.lme2
#' @rdname sCorrect
#' @export
sCorrect.lme2 <- function(object, ...){
    class(object) <- setdiff(class(object),"lme2")
    return(sCorrect(object, ...))    
}

## * sCorrect.lvmfit
#' @rdname sCorrect
#' @export
sCorrect.lvmfit <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                            score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                            param = NULL, data = NULL,
                            tol = 1e-5, n.iter = 20, trace = 0,
                            ...){

### ** Check valid lvm object
    if("multigroupfit" %in% class(object)){
        stop("sCorrect does not handle multigroup models \n")
    }
    ## if(!is.null(object$cluster)){
    ##     stop("sCorrect does not handle lvmfit object with cluster \n")
    ## }
   
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
    
### ** Extract quantities from object
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

    name.latent <- latent(object)
    n.latent <- length(name.latent)

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
        cat("* Compute conditional moments and their derivative ")
    }
    object$conditionalMoment <- conditionalMoment(object, data = data, param = model.param,
                                                  first.order = TRUE, second.order = FALSE, usefit = TRUE)
    if(df == TRUE && (numeric.derivative == FALSE)){
        object$conditionalMoment$d2Moment.init <- skeletonDtheta2(lava::Model(object),
                                                                  data = data,
                                                                  df.param.all = object$conditionalMoment$df.param,
                                                                  param2originalLink = object$conditionalMoment$param2originalLink,
                                                                  name.latent = name.latent)
    }
    if(trace>0){
        cat("- done \n")
    }

#### ** Compute residuals
    if(trace>0){
        cat("* Extract residuals ")
    }
    epsilon <- as.matrix(data[, name.endogenous,drop=FALSE] - object$conditionalMoment$mu)
    ## residuals(object) - epsilon
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

### ** args
    args <- list(adjust.Omega = adjust.Omega,
                 adjust.n = adjust.n,                     
                 df = df,
                 numeric.derivative = numeric.derivative,
                 tol = tol, n.iter = n.iter)
    if(df && numeric.derivative){
        argsNumDeriv <- list(data = data)
    }else{
        argsNumDeriv <- list()
    }
    
### ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = epsilon,
                     name.param = name.param,
                     name.endogenous = name.endogenous,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     score = score,
                     derivative = derivative,
                     args = args,
                     argsNumDeriv = argsNumDeriv,
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
.sCorrect <- function(object, param, epsilon, 
                      name.param, name.endogenous, 
                      n.cluster, index.Omega,
                      score, derivative, args, argsNumDeriv, trace){

    n.param <- length(param)
    if(!is.null(index.Omega)){
        n.endogenous.cluster <- lapply(index.Omega,length)        
    }else{
        n.endogenous.cluster <- NULL
    }
    name.3deriv <- object$conditionalMoment$name.3deriv
    
    ## ** check names of the mean and variance parameters
    name.meanparam <- names(object$conditionalMoment$dmu)
    name.meanparam <- as.character(sort(factor(name.meanparam, levels = name.param)))
    if(any(is.na(name.meanparam))){
        stop("An element in name.meanparam is not in name.param. \n")
    }
    if(length(name.meanparam)>0 && !identical(sort(name.meanparam),sort(names(object$conditionalMoment$dmu)))){
        stop("Mismatch first derivative of the conditional mean and name.meanparam \n")
    }

    name.varparam <- names(object$conditionalMoment$dOmega)
    name.varparam <- as.character(sort(factor(name.varparam, levels = name.param)))
    if(any(is.na(name.varparam))){
        stop("An element in name.varparam is not in name.param. \n")
    }
    if(length(name.varparam)>0 && !identical(sort(name.varparam),sort(names(object$conditionalMoment$dOmega)))){
        stop("Mismatch first derivative of the conditional variance and name.varparam \n")
    }
    if(length(name.varparam)==0){
        args$adjust.n <- FALSE
        args$adjust.Omega <- FALSE
    }
    
    ## ** corrected ML estimates
    object  <- .estimate2(object = object,
                          epsilon = epsilon,
                          n.cluster = n.cluster,
                          name.param = name.param,
                          name.meanparam = name.meanparam,
                          name.varparam = name.varparam,
                          name.endogenous = name.endogenous,
                          index.Omega = index.Omega, ## mode2
                          adjust.Omega = args$adjust.Omega,
                          adjust.n = args$adjust.n,
                          tol = args$tol, n.iter = args$n.iter,
                          trace = trace)

    Omega <- object$conditionalMoment$Omega    
    if(!is.null(index.Omega)){
        OmegaM1 <- lapply(1:n.cluster, function(iC){
            return(solve(Omega[index.Omega[[iC]],index.Omega[[iC]]]))
        })    
    }else{
        OmegaM1 <- chol2inv(chol(Omega))
    }
    
    ## ** corrected score
    if(score){
        if(trace>0){
            if(args$adjust.n == FALSE && args$adjust.Omega == FALSE){
                cat("* Compute score ")
            }else{
                cat("* Compute corrected score ")
            }
        }
        object$dVcov$score <- .score2(epsilon = object$dVcov$residuals,
                                      Omega = Omega,
                                      OmegaM1 = OmegaM1,
                                      dmu = object$conditionalMoment$dmu,
                                      dOmega = object$conditionalMoment$dOmega,
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
    if(args$df == FALSE || length(name.3deriv)==0){
        object$dVcov$dVcov.param <- NULL
    }else if(derivative == "none"){
        object$dVcov$dVcov.param <- NA
    }else if(derivative == "numeric"){
        if(trace>0){
            cat("Compute first derivative of the information matrix using numerical differentiation ")
        }
        if(args$adjust.Omega || args$adjust.n){
            warning("The numerical derivative of the information matrix is computed ignoring the small sample correction \n")
        }

        args.tempo <- args
        args.tempo$data <- argsNumDeriv$data
        args.tempo$df <- FALSE
        args.tempo$score <- FALSE

        ## *** direct computation of the variance-covariance matrix
        calcVcov <- function(iParam){ # x <- p.obj
            pp <- param
            pp[names(iParam)] <- iParam
            vcov.param <- do.call(sCorrect,
                                  args = c(list(object, param = pp), args.tempo))$vcov.param
            return(vcov.param)
            ## return(solve(vcov.param))
        }

        ### *** numerical derivative
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }
        jac.param <- param[name.3deriv]
        res.numDeriv <- numDeriv::jacobian(calcVcov, x = jac.param, method = "Richardson")
        
        object$dVcov$dVcov.param <- array(res.numDeriv,
                                          dim = c(n.param,n.param,length(name.3deriv)),
                                          dimnames = list(name.param, name.param, name.3deriv))
        if(trace>0){
            cat("- done \n")
        }
        
    }else if(derivative == "analytic"){
        if(trace>0){
            cat("* Compute first derivative of the information matrix using analytic formula ")
        }

        ## update conditional moments
        resD2 <- skeletonDtheta2(object)

        dInfo.dtheta <- .dInformation2(dmu = object$conditionalMoment$dmu,
                                       d2mu = resD2$d2mu,
                                       dOmega = object$conditionalMoment$dOmega,
                                       d2Omega = resD2$d2Omega,
                                       Omega = Omega,
                                       OmegaM1 = OmegaM1,
                                       n.corrected = object$dVcov$n.corrected,
                                       n.cluster = n.cluster,
                                       index.Omega = index.Omega,
                                       leverage = object$dVcov$leverage,
                                       name.param  = name.param,
                                       name.3deriv = name.3deriv)
        p3 <- dim(dInfo.dtheta)[3]
        object$dVcov$dVcov.param <- array(NA, dim = dim(dInfo.dtheta), dimnames = dimnames(dInfo.dtheta))
        
        for(iP in 1:p3){
            object$dVcov$dVcov.param[,,iP] <- - object$dVcov$vcov.param %*% dInfo.dtheta[,,iP] %*% object$dVcov$vcov.param
            ## object$dVcov$dVcov.param[,,iP] <- dInfo.dtheta[,,iP]
        }

        if(trace>0){
            cat("- done \n")
        }
    }
       
    ## ** export
    object$dVcov$args <- args
    return(object$dVcov)
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
## * sCorrect<-.lm2
#' @rdname sCorrect
#' @export
`sCorrect<-.lm2` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
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
## * sCorrect<-.gls2
#' @rdname sCorrect
#' @export
`sCorrect<-.gls2` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
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
## * sCorrect<-.lme2
#' @rdname sCorrect
#' @export
`sCorrect<-.lme2` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    return(x)
}    
## * sCorrect<-.lvmfit
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit` <- function(x, ..., value){
    dots <- list(...)
    safeMode <- dots$safeMode
    dots[["safeMode"]] <- NULL
    
    if(identical(safeMode,TRUE)){
        x$sCorrect <-  try(do.call(sCorrect,
                                   args = c(list(x, adjust.Omega = value, adjust.n = value),
                                            dots) ), silent = TRUE)
        if(value == TRUE && inherits(x$sCorrect,"try-error")){
            warn <- x$sCorrect
            x$sCorrect <- do.call(sCorrect,
                                  args = c(list(x, adjust.Omega = value, adjust.n = FALSE),
                                           dots) )
            attr(x$sCorrect,"warning") <- warn
            warning("sCorrect failed and has been re-run setting the argument \'adjust.n\' to FALSE \n",
                    "see the attribute \"warning\" of object$sCorrect for the error message \n")
        }
    }else{
        x$sCorrect <-  do.call(sCorrect,
                               args = c(list(x, adjust.Omega = value, adjust.n = value),
                                        dots) )
    }
    class(x) <- append("lvmfit2",class(x))

    return(x)
}    

## * sCorrect<-.lvmfit2
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit2` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    return(x)
}



##----------------------------------------------------------------------
### sCorrect.R ends here








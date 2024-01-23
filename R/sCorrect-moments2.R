### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: jan 23 2024 (10:28) 
##           By: Brice Ozenne
##     Update #: 1955
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * moments2 - documentation
#' @title Compute Key Quantities of a Latent Variable Model
#' @description Compute conditional mean, conditional variance, their first and second derivative regarding model parameters, as well as various derivatives of the log-likelihood.
#' @name moments2
#' 
#' @param object a latent variable model.
#' @param data [data.frame] dataset if different from the one used to fit the model.
#' @param param [numeric vector] value of the model parameters if different from the estimated ones.
#' @param initialize [logical] Pre-compute quantities dependent on the data but not on the parameters values.
#' @param usefit [logical] Compute key quantities based on the parameter values.
#' @param update.dmoment [logical] should the first derivative of the moments be computed/updated?
#' @param update.d2moment [logical] should the second derivative of the the moments be computed/updated?
#' @param score [logical] should the score be output?
#' @param information [logical] should the expected information be output?
#' @param hessian [logical] should the hessian be output?
#' @param vcov [logical] should the variance-covariance matrix based on the expected information be output?
#' @param dVcov [logical] should the derivative of the variance-covariance matrix be output?
#' @param dVcov.robust [logical]  should the derivative of the robust variance-covariance matrix be output?
#' @param Psi [matrix]  Average first order bias in the residual variance. Only necessary for computing adjusted residuals.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient (\code{conditionalMoment.lvm}).
#' \item an advanced one that require the model coefficients (\code{conditionalMoment.lvmfit}). 
#' }
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' d <- lava::sim(m,1e2)
#' e <- estimate(m, d)
#'
#' ## basic pre-computation
#' res1 <- moments2(e, data = d, initialize = TRUE, usefit = FALSE,
#'                 score = TRUE, information = TRUE, hessian = TRUE, vcov = TRUE,
#'                 dVcov = TRUE, dVcov.robust = TRUE, residuals = TRUE, leverage = FALSE,
#'                 derivative = "analytic")
#' res1$skeleton$param$Sigma
#' 
#' ## full pre-computation
#' res2 <- moments2(e, param = coef(e), data = d, initialize = TRUE, usefit = TRUE,
#'                 score = TRUE, information = TRUE, hessian = TRUE, vcov = TRUE,
#'                 dVcov = TRUE, dVcov.robust = TRUE, residuals = TRUE, leverage = FALSE,
#'                 derivative = "analytic")
#' res2$moment$Omega
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
#' @export
`moments2` <-
    function(object, param, data, weights, Omega, Psi,
             initialize, usefit,
             update.dmoment, update.d2moment, score, information, hessian, vcov, dVcov, dVcov.robust, residuals, leverage, derivative) UseMethod("moments2")

## * moments2.lvm
#' @rdname moments2
#' @export
moments2.lvm <- function(object, param = NULL, data = NULL, weights = NULL, Omega = NULL, Psi = NULL,
                         initialize = TRUE, usefit = TRUE,
                         update.dmoment = TRUE, update.d2moment = TRUE, score = TRUE, information = TRUE, hessian = TRUE, vcov = TRUE, dVcov = TRUE, dVcov.robust = TRUE, residuals = TRUE, leverage = TRUE,
                         derivative = "analytic"){
    if(lava.options()$debug){cat("moments2 \n")}
    
    ## ** sanity checks
    if(initialize == FALSE && is.null(object$sCorrect)){
        stop("Initialization of the moments missing. \n",
             "Consider setting the argument \'initialize\' to TRUE \n")
    }
    if(!missing(derivative)){
        derivative <- match.arg(derivative, choices = c("analytic","numeric"))
    }
    
    ## ** initialize
    if(score || information || hessian || vcov || dVcov || dVcov.robust){
        first.order <- TRUE
    }else{
        first.order <- FALSE
    }
    if(hessian || dVcov || dVcov.robust){
        second.order <- TRUE
    }else{
        second.order <- FALSE
    }

    ## NOTE: the estimation of the leverage depends on the information matrix and vice versa.
    ##       moments2 is a single run function called by estimate2 which take care of iterating until reaching a stable point
    previous.vcov.param <- object$sCorrect$vcov.param

    ## ** skeleton
    if(initialize){
        out <- list()

        ## *** get information from object
        out$endogenous <- lava::endogenous(object, format = "wide")
        out$latent <- lava::latent(object)

        if(is.null(data)){
            out$data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                                    envir = parent.env(environment()), rm.na = TRUE)
        }else{
            out$data <- as.data.frame(data)
        }
        out$X <- as.matrix(out$data[,lava::manifest(object),drop=FALSE])
        out$cluster <- .getGroups2(object, data = out$data, endogenous = out$endogenous)

        reserved.name <- c("XXvalueXX","XXendogenousXX","XXendogenousXX.index","XXclusterXX")
        if(any(colnames(out$X) %in% reserved.name)){
            stop("\"",paste(reserved.name[colnames(out$X) %in% reserved.name], collapse="\" \""),"\" should not correspond to variable names \n",
                 "It is used internally for data manipulation \n")
        }

        ## *** reshape dataset (convert to long format)
        X.latent <- matrix(NA, nrow = NROW(out$X), ncol = length(out$latent),
                           dimnames = list(NULL, out$latent))

        X.long <- reshape2::melt(data.frame(XXclusterXX = out$cluster$name.cluster, out$X, X.latent),
                                 id.vars = "XXclusterXX",
                                 measure.vars = c(out$endogenous,out$latent),
                                 variable.name = "XXendogenousXX",
                                 value.name = "XXvalueXX")
        X.long <- cbind(X.long,out$X[X.long$XXclusterXX,,drop=FALSE]) ## add all variable (endo and exo) in case some endo are regressors

        X.wide <- out$X[,out$endogenous,drop=FALSE]

        out$old2new.order <- X.long[,c("XXclusterXX","XXendogenousXX","XXvalueXX")]
        names(out$old2new.order) <- c("XXclusterXX.old","XXendogenousXX.old","XXvalueXX.old")
        X.long <- X.long[order(X.long$XXclusterXX,X.long$XXendogenousXX),]
        out$old2new.order$XXclusterXX.new <- out$old2new.order$XXclusterXX
        out$old2new.order$XXendogenousXX.new <- out$old2new.order$XXendogenousXX
        out$old2new.order$XXvalueXX.new <- out$old2new.order$XXvalueXX
                
        ## *** identify missing pattern
        pattern <- X.wide
        pattern[!is.na(pattern)] <- 1
        pattern[is.na(pattern)] <- 0

        if(!is.null(object$call$missing) && object$call$missing==FALSE){
            unique.pattern <- unique(pattern[rowSums(pattern==0)==0,,drop=FALSE])
            if(length(unique.pattern)==0){
                stop("All clusters contain at least one missing value - cannot perform complete case analysis. \n",
                     "Consider setting the argument \'missing\' to TRUE. \n")
            }
            XXclusterXX.NA <- unique(X.long$XXclusterXX[rowSums(is.na(X.long[,out$endogenous,drop=FALSE]))>0])
            X.long[X.long$XXclusterXX %in% XXclusterXX.NA,c("XXvalueXX",out$endogenous)] <- NA ## add missing value to the other observation of the same cluster
        }else{
            unique.pattern <- unique(pattern[rowSums(pattern==1)>0,,drop=FALSE])
        }
        
        name.pattern <- apply(unique.pattern, MARGIN = 1, FUN = paste0, collapse = "")
        rownames(unique.pattern) <- name.pattern

        out$missing$pattern <- tapply(1:out$cluster$n.cluster,apply(pattern, MARGIN = 1, FUN = paste0, collapse=""),list)[name.pattern]
        out$missing$unique.pattern <- unique.pattern
        out$missing$name.pattern <- name.pattern

        ## *** initialize conditional moments
        out$skeleton <- skeleton(object, X = X.long,
                                 endogenous = out$endogenous, latent = out$latent,
                                 n.cluster = out$cluster$n.cluster,
                                 index.Omega = out$cluster$index.Omega)

        ## *** initialize partial derivatives of the conditional moments
        out$skeleton <- skeletonDtheta(out$skeleton,
                                       X = X.long,
                                       endogenous = out$endogenous, latent = out$latent,
                                       missing.pattern = out$missing$pattern,
                                       unique.pattern = out$missing$unique.pattern,
                                       name.pattern = out$missing$name.pattern,
                                       n.cluster = out$cluster$n.cluster,
                                       index.Omega = out$cluster$index.Omega)

        ## *** initialize second order partial derivatives of the conditional moments
        ## GS <- skeletonDtheta2(out$skeleton)
        out$skeleton <- skeletonDtheta2(out$skeleton)

        ## *** weights
        if(!is.null(object$weights)){
            out$weights <- object$weights[,1]
        }

    }else{
        ## subset to remove existing results
        rm.name <- c("moment","dmoment","d2moment","score","vcov.param","information","hessian","dInformation","dVcov.param","dRvcov.param","leverage","residuals")
        out <- object$sCorrect[setdiff(names(object$sCorrect),rm.name)]
    }

    ## ** update according to the value of the model coefficients
    if(usefit){
        if(is.null(param)){
            param.tempo <- stats::coef(object, type = 2, labels = 1)
            out$param <- stats::setNames(param.tempo[,"Estimate"],rownames(param.tempo))[out$skeleton$Uparam]
            ## out$name.param <- out$name.param[names(stats::coef(object))]
        }else{
            if(all(names(param) %in% out$skeleton$Uparam)){ ## using user-defined names
                ## e.g. mu1 mu2 Y1~X1 Y2~X1 sigma
                out$param[names(param)] <- param
            }else if(all(names(param) %in% names(out$skeleton$originalLink2param))){ ## using original link
                ## e.g. Y1 Y2 Y1~X1 Y2~X1 Y1~~Y1
                out$param[out$skeleton$Uparam[match(names(param),names(out$skeleton$originalLink2param))]] <- param
            }else{
                stop("Could not find model parameter(s) corresponding to the name(s): \"",paste(setdiff(names(param),c(out$skeleton$Uparam,names(out$skeleton$Uparam))), collapse="\" \""),"\" \n")
            }
        }
        if(is.null(weights)){
            weights <- object$weights[,1]
        }
        
        out$name.param <- stats::setNames(out$skeleton$type[!is.na(out$skeleton$type$lava),"param"],
                                   out$skeleton$type[!is.na(out$skeleton$type$lava),"originalLink"])

        ## *** conditional moments
        out$moment <- updateMoment(skeleton = out$skeleton$param,
                                   value = out$skeleton$value,
                                   toUpdate = out$skeleton$toUpdate.moment,
                                   param = out$param, Omega = Omega,
                                   name.pattern = out$missing$name.pattern,
                                   unique.pattern = out$missing$unique.pattern,
                                   endogenous = out$endogenous,
                                   latent = out$latent,
                                   n.cluster = out$cluster$n.cluster)

        ## *** first order derivatives
        if(update.dmoment && first.order){
            out$dmoment <- updateDMoment(moment = out$moment,
                                         skeleton = out$skeleton,
                                         param = out$param)
        }

        
        ## *** second order derivatives
        if(update.d2moment && second.order){
            out$d2moment <- updateD2Moment(moment = out$moment,
                                           skeleton = out$skeleton,
                                           param = out$param)
        }
        

        ## *** update residuals
        if(residuals || score || information || vcov || hessian || leverage){
            if(!is.null(attr(Omega,"Omega.residuals"))){
                OOmega <- attr(Omega,"Omega.residuals")
            }else{
                OOmega <- out$moment$Omega
            }

            out$residuals <- .adjustResiduals(epsilon = out$skeleton$param$endogenous - out$moment$mu,
                                              Omega = OOmega, Psi = Psi, ## Note: if Psi is null returns epsilon i.e. no adjustment
                                              name.pattern = out$missing$name.pattern, missing.pattern = out$missing$pattern, unique.pattern = out$missing$unique.pattern,
                                              endogenous = out$endogenous, n.cluster = out$cluster$n.cluster)

        }
        ## mean(out$residuals^2)
        ## out$moment$Omega
        
        ## *** score
        if(score){
            out$score <- .score2(dmu = out$dmoment$dmu,
                                 dOmega = out$dmoment$dOmega,                    
                                 epsilon = out$residuals,
                                 OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                 missing.pattern = out$missing$pattern,
                                 unique.pattern = out$missing$unique.pattern,
                                 name.pattern = out$missing$name.pattern,
                                 name.param = out$skeleton$Uparam,
                                 name.meanparam = out$skeleton$Uparam.mean,
                                 name.varparam = out$skeleton$Uparam.var,
                                 n.cluster = out$cluster$n.cluster,
                                 weights = weights)
        }

        ## *** leverage
        if((leverage || information || vcov) && !is.null(previous.vcov.param)){
            if(!is.null(attr(Omega,"Omega.leverage"))){
                OOmega <- attr(Omega,"Omega.leverage")
            }else{
                OOmega <- out$moment$Omega
            }
            if(!is.null(attr(Omega,"dOmega.leverage"))){
                ddOOmega <- attr(Omega,"dOmega.leverage")
            }else{
                ddOOmega <- out$dmoment$dOmega
            }

            out$leverage <- .leverage2(Omega = OOmega, 
                                       epsilon = out$residuals,
                                       dmu = aperm(abind::abind(out$dmoment$dmu, along = 3), perm = c(3,2,1)),
                                       dOmega = ddOOmega,
                                       vcov.param = previous.vcov.param,
                                       name.pattern = out$missing$name.pattern,
                                       missing.pattern = out$missing$pattern,
                                       unique.pattern = out$missing$unique.pattern,
                                       endogenous = out$endogenous,
                                       n.endogenous = length(out$endogenous),
                                       param = out$skeleton$Uparam,
                                       param.mean = out$skeleton$Uparam.mean,
                                       param.var = out$skeleton$Uparam.var,
                                       n.cluster = out$cluster$n.cluster)
        }

        ## *** information matrix
        if(information || vcov){
            if(!is.null(attr(Omega,"dOmega.leverage"))){
                ddOOmega <- attr(Omega,"dOmega.leverage")
            }else{
                ddOOmega <- out$dmoment$dOmega
            }

            out$information <- .information2(dmu = out$dmoment$dmu,
                                             dOmega = ddOOmega,##out$dmoment$dOmega,
                                             OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                             missing.pattern = out$missing$pattern,
                                             unique.pattern = out$missing$unique.pattern,
                                             name.pattern = out$missing$name.pattern,
                                             grid.mean = out$skeleton$grid.dmoment$mean, 
                                             grid.var = out$skeleton$grid.dmoment$var, 
                                             name.param = out$skeleton$Uparam,
                                             leverage = out$leverage,
                                             n.cluster = out$cluster$n.cluster,
                                             weights = weights)
        }

        if(vcov){
            out$vcov.param  <- .info2vcov(out$information, attr.info = FALSE)
        }
    
        ## *** hessian
        if(hessian || (derivative == "analytic") && dVcov.robust){
            out$hessian <- .hessian2(dmu = out$dmoment$dmu,
                                     dOmega = out$dmoment$dOmega,
                                     d2mu = out$d2moment$d2mu,
                                     d2Omega = out$d2moment$d2Omega,
                                     epsilon = out$residuals,                                     
                                     OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                     missing.pattern = out$missing$pattern,
                                     unique.pattern = out$missing$unique.pattern,
                                     name.pattern = out$missing$name.pattern,
                                     grid.mean = out$skeleton$grid.dmoment$mean, 
                                     grid.var = out$skeleton$grid.dmoment$var, 
                                     grid.hybrid = out$skeleton$grid.dmoment$hybrid, 
                                     name.param = out$skeleton$Uparam,
                                     leverage = out$leverage,
                                     n.cluster = out$cluster$n.cluster,
                                     weights = weights)
        }

        ## *** dVcov.param (model based variance, analytic)
        if((dVcov || dVcov.robust) && (derivative == "analytic")){
            out$dInformation <- .dInformation2(dmu = out$dmoment$dmu,
                                               dOmega = out$dmoment$dOmega,
                                               d2mu = out$d2moment$d2mu,
                                               d2Omega = out$d2moment$d2Omega,
                                               OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                               missing.pattern = out$missing$pattern,
                                               unique.pattern = out$missing$unique.pattern,
                                               name.pattern = out$missing$name.pattern,
                                               grid.3varD1 = out$skeleton$grid.3varD1,
                                               grid.2meanD1.1varD1 = out$skeleton$grid.2meanD1.1varD1,
                                               grid.2meanD2.1meanD1 = out$skeleton$grid.2meanD2.1meanD1,
                                               grid.2varD2.1varD1 = out$skeleton$grid.2varD2.1varD1,
                                               name.param = out$skeleton$Uparam,
                                               leverage = out$leverage,
                                               n.cluster = out$cluster$n.cluster,
                                               weights = weights)

            ## delta method
            out$dVcov.param <- .dVcov.param(vcov.param = out$vcov.param,
                                            dInformation = out$dInformation,
                                            n.param = length(out$skeleton$Uparam),
                                            name.param = out$skeleton$Uparam)
        }

        ## *** dRvcov.param  (robust variance, analytic)
        if(dVcov.robust && (derivative == "analytic")){

        out$dRvcov.param <- .dRvcov.param(score = out$score,
                                          hessian = out$hessian,
                                          vcov.param = out$vcov.param,
                                          dVcov.param = out$dVcov.param,
                                          n.param = length(out$skeleton$Uparam),
                                          name.param = out$skeleton$Uparam)

    }

    ## *** dVcov.param and dRvcov.param (numeric derivatives)
    if((dVcov || dVcov.robust) && derivative == "numeric"){
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }

        ## range(out$score - .warper.numDev(value = out$param, object = object, type = "score"))
        ## range(out$hessian - .warper.numDev(value = out$param, object = object, type = "hessian"))
        ## range(out$vcov.param - .warper.numDev(value = out$param, object = object, type = "vcov.model"))

        param <- out$param
        name.param <- names(out$param)
        n.param <- length(param)
        n.cluster <- out$cluster$n.cluster
        object2 <- object
        object2$sCorrect <- out

        ## *** hessian
        ## print(range(.warper.numDev(param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "score")-out$score))
        num.hessian <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "score", method = "Richardson")

        out$hessian <- aperm(array(num.hessian, dim = c(n.cluster,n.param,n.param),
                                               dimnames = list(NULL, name.param, name.param)), perm = 3:1)
        
        ## *** dInformation
        ## print(range(.warper.numDev(param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "information")-out$information))
        num.information <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "information", method = "Richardson")

        out$dInformation <- array(num.information, dim = c(n.param,n.param,n.param),
                                              dimnames = list(name.param, name.param, name.param))

        ## *** dVcov.param
        ## print(range(.warper.numDev(param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "vcov")-out$vcov.param))
        num.dVcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "vcov", method = "Richardson")

        out$dVcov.param <- array(num.dVcov.param, dim = c(n.param,n.param,n.param),
                                             dimnames = list(name.param, name.param, name.param))


        ## *** dRvcov.param
        num.dRvcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, weights = weights, Omega = Omega, Psi = Psi, type = "vcov.robust", method = "Richardson")

        out$dRvcov.param <- array(num.dRvcov.param, dim = c(n.param,n.param,n.param),
                                              dimnames = list(name.param, name.param, name.param))
    }        
    
        
    }

    ## ** export
    return(out)
}

## * moments2.lvmfit
#' @rdname moments2
#' @export
moments2.lvmfit <- moments2.lvm

## * .wraper.numDev (helper)
.warper.numDev <- function(value, object, type, weights, Psi = NULL, Omega = NULL){ # x <- p.obj
    ## CANNOT DO DIRECTLY VIA moments2
    ## (because Omega should not be fixed in updateMoment but it should be fixed (as well as Psi) when updating the residuals)
    
    out <- object$sCorrect
    out$moment <- updateMoment(skeleton = out$skeleton$param,
                               value = out$skeleton$value,
                               toUpdate = out$skeleton$toUpdate.moment,
                               param = value, Omega = NULL,
                               name.pattern = out$missing$name.pattern,
                               unique.pattern = out$missing$unique.pattern,
                               endogenous = out$endogenous,
                               latent = out$latent,
                               n.cluster = out$cluster$n.cluster)

    out$dmoment <- updateDMoment(moment = out$moment,
                                 skeleton = out$skeleton,
                                 param = out$param)

    if(type %in% c("score","vcov.robust")){
        out$residuals <- .adjustResiduals(epsilon = out$skeleton$param$endogenous - out$moment$mu,
                                          Omega = Omega, Psi = Psi, ## FIXED
                                          name.pattern = out$missing$name.pattern, missing.pattern = out$missing$pattern, unique.pattern = out$missing$unique.pattern,
                                          endogenous = out$endogenous, n.cluster = out$cluster$n.cluster)

        out$score <- .score2(dmu = out$dmoment$dmu,
                             dOmega = out$dmoment$dOmega,                    
                             epsilon = out$residuals,
                             OmegaM1 = out$moment$OmegaM1.missing.pattern,
                             missing.pattern = out$missing$pattern,
                             unique.pattern = out$missing$unique.pattern,
                             name.pattern = out$missing$name.pattern,
                             name.param = out$skeleton$Uparam,
                             name.meanparam = out$skeleton$Uparam.mean,
                             name.varparam = out$skeleton$Uparam.var,
                             n.cluster = out$cluster$n.cluster,
                             weights = weights)
    }

    if(type %in% c("information","vcov","vcov.robust")){
        out$information <- .information2(dmu = out$dmoment$dmu,
                                         dOmega = out$dmoment$dOmega,
                                         OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                         missing.pattern = out$missing$pattern,
                                         unique.pattern = out$missing$unique.pattern,
                                         name.pattern = out$missing$name.pattern,
                                         grid.mean = out$skeleton$grid.dmoment$mean, 
                                         grid.var = out$skeleton$grid.dmoment$var, 
                                         name.param = out$skeleton$Uparam,
                                         leverage = out$leverage, ## FIXED
                                         n.cluster = out$cluster$n.cluster,
                                         weights = weights)
    }

    if(type %in% c("vcov","vcov.robust")){
        out$vcov  <- .info2vcov(out$information, attr.info = FALSE)
    }

    if(type %in% c("vcov.robust")){
        out$vcov.robust  <- out$vcov %*% crossprod(out$score) %*% out$vcov
    }


    return(out[[type]])
}


## * .info2vcov (helper)
#' @title Inverse the Information Matrix
#' @description Compute the inverse of the information matrix.
#' @name vcov2-internal
#'
#' @param information [matrix] information matrix to be inverted.
#' @param attr.info [logical] should the information matrix be returned as an attribute?
#' 
#' @keywords internal
.info2vcov <- function(information, attr.info = FALSE){
    vcov <- try(chol2inv(chol(information)), silent = TRUE)
    if(inherits(vcov, "try-error")){
        vcov <- try(solve(information), silent = TRUE)
        if(inherits(vcov, "try-error")){ ## try by block
            cat("Singular information matrix: try to inverse it by block \n")
            information.N0 <- abs(information)>1e-10
            remaining.var <- colnames(information)
            vcov <- matrix(0, nrow = NROW(information), ncol = NCOL(information),
                           dimnames = dimnames(information))
            while(length(remaining.var)>0){
                current.set <- remaining.var[1]
                new.set <- unique(unlist(apply(information.N0[current.set,,drop=FALSE],1,function(iRow){list(names(iRow[iRow==1]))})))
                while(length(current.set)<length(new.set)){
                    current.set <- new.set
                    new.set <- unique(unlist(apply(information.N0[current.set,,drop=FALSE],1,function(iRow){list(names(iRow[iRow==1]))})))
                }
                if(length(new.set)>0){
                    iTry <- try(solve(information[current.set,current.set]), silent = TRUE)
                    if(inherits(iTry,"try-error")){
                        vcov[current.set,current.set] <- NA
                    }else{
                        vcov[current.set,current.set] <- iTry
                    }
                }
                remaining.var <- setdiff(remaining.var, current.set)
            }
        }
    }
    if(attr.info){
        attr(vcov,"information") <- information
    }
    if(!inherits(vcov, "try-error")){
        dimnames(vcov) <- dimnames(information)
    }
    return(vcov)
}

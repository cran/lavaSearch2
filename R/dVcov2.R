### dVcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: jan 19 2018 (15:53) 
##           By: Brice Ozenne
##     Update #: 225
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - dVcov2
#' @title  Compute the Derivative of the Information Matrix
#' @description Compute the derivative of the information matrix.
#' @name dVcov2
#'
#' @param object a gls, lme, or lvm object.
#' @param x same as object.
#' @param cluster the grouping variable relative to which the observations are iid.
#'                Only required for gls models with no correlation argument.
#' @param vcov.param the variance-covariance matrix of the estimates.
#' @param adjust.residuals Small sample correction: should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.
#' @param value same as adjust.residuals.
#' @param numericDerivative If TRUE, the degree of freedom are computed using a numerical derivative.
#' @param return.score [for internal use] export the score.
#' @param ... arguments to be passed to lower level methods.
#' 
#' @export
`dVcov2` <-
  function(object, ...) UseMethod("dVcov2")


## * dVcov.lm
#' @rdname dVcov2
#' @export
dVcov2.lm <- function(object, adjust.residuals = FALSE,
                      return.score = FALSE, ...){

    ## ** extract information
    X <- stats::model.matrix(object)
        
    epsilonD2 <- residuals2(object, adjust.residuals = adjust.residuals,
                            return.vcov.param = TRUE)
    p <- c(stats::coef(object), sigma2 = mean(stats::residuals(object)^2))
    name.param <- names(p)
    n.param <- length(p)
    vcov.param <- attr(epsilonD2, "vcov.param")
    sigma2.cor <- attr(epsilonD2, "sigma2")

    dVcov.dtheta <- array(0,dim = c(n.param, n.param, 1),
                          dimnames = list(name.param, name.param, "sigma2"))

    dVcov.dtheta[setdiff(name.param,"sigma2"),setdiff(name.param,"sigma2"),"sigma2"] <- solve(t(X) %*% X)
    dVcov.dtheta["sigma2","sigma2","sigma2"] <- 2*vcov.param["sigma2","sigma2"]/sigma2.cor

    ## ** score
    if(return.score){
        attr(epsilonD2, "vcov.param") <- NULL
        attr(epsilonD2, "sigma2") <- NULL
        X <- stats::model.matrix(object)
        out.score <- cbind(sweep(X, MARGIN = 1, FUN = "*",  STATS = epsilonD2) / sigma2.cor,
                           -1/(2*sigma2.cor) + 1/(2*sigma2.cor^2) * epsilonD2^2)

        attr(dVcov.dtheta, "score") <- out.score
    }
    
    
    ## ** export
    attr(dVcov.dtheta, "vcov.param") <- vcov.param
    attr(dVcov.dtheta, "param") <- p
    return(dVcov.dtheta)
    
}
## * dVcov2.gls
#' @rdname dVcov2
#' @export
dVcov2.gls <- function(object, cluster, vcov.param = NULL,
                       adjust.residuals = FALSE, numericDerivative = FALSE,
                       return.score = FALSE, ...){

    p <- .coef2(object)
    data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE)

    n.param <- length(p)
    name.param <- names(p)

    as.clubSandwich <- 1

    ## ** param with non-zero third derivative
    keep.param <- setdiff(name.param, attr(.coef2(object),"mean.coef"))

    ## ** pre-compute quantities
    if(return.score || numericDerivative == FALSE){
        epsilonD2 <- residuals2(object, p = p, data = data, cluster = cluster,
                                adjust.residuals = adjust.residuals,
                                as.clubSandwich = as.clubSandwich,
                                return.prepareScore2 = TRUE, return.vcov.param = TRUE, second.order = TRUE)
        pS2  <- attr(epsilonD2, "prepareScore2")
    }
    
    ### ** Compute the gradient
    if(numericDerivative){
        
        ### *** Define function to compute the variance-covariance matrix
        calcSigma <- function(iParam){ # x <- p.obj
            pp <- p
            pp[names(iParam)] <- iParam
            res.tempo <- residuals2(object, cluster = cluster, p = pp, data = data,
                                    adjust.residuals = adjust.residuals,
                                    as.clubSandwich = as.clubSandwich,
                                    return.vcov.param = TRUE, second.order = FALSE)
            vcov.tempo <- attr(res.tempo, "vcov.param")
            attr(vcov.param, "warning")  <- attr(res.tempo, "warning")
            return(vcov.tempo)
        }
        
        if(is.null(vcov.param)){
            vcov.param <- calcSigma(p)
        }

        ### *** numerical derivative
        jac.param <- p[keep.param]
        res.tempo <- numDeriv::jacobian(calcSigma, x = jac.param, method = "Richardson")

        dVcov.dtheta <- array(res.tempo,
                              dim = c(n.param,n.param,length(jac.param)),
                              dimnames = list(name.param, name.param, keep.param))

    }else{        
        vcov.param  <- attr(epsilonD2, "vcov.param")
        attr(vcov.param, "warning")  <- attr(epsilonD2, "warning")
        
        dInfo.dtheta <- .dinformation2(dmu.dtheta = pS2$dmu.dtheta,
                                       d2mu.dtheta2 = NULL,
                                       dOmega.dtheta = pS2$dOmega.dtheta,
                                       d2Omega.dtheta2 = pS2$d2Omega.dtheta2,
                                       Omega = pS2$Omega,
                                       ls.indexOmega = pS2$ls.indexOmega,
                                       hat = pS2$hat,
                                       n.param  = pS2$n.param,
                                       name.param  = pS2$name.param,
                                       name.deriv = keep.param,
                                       n.cluster = pS2$n.cluster)

        p3 <- dim(dInfo.dtheta)[3]
        dVcov.dtheta <- array(NA, dim = dim(dInfo.dtheta), dimnames = dimnames(dInfo.dtheta))
        for(iP in 1:p3){
            dVcov.dtheta[,,iP] <- - vcov.param %*% dInfo.dtheta[,,iP] %*% vcov.param
        }
        
    }

    ## ** score
    if(return.score){
        attr(epsilonD2, "prepareScore2") <- NULL
        attr(epsilonD2, "vcov.param") <- NULL
        
        out.score <- .score2(dmu.dtheta = pS2$dmu.dtheta,
                             dOmega.dtheta = pS2$dOmega.dtheta,
                             epsilon = epsilonD2,
                             Omega = pS2$Omega,                         
                             ls.indexOmega = pS2$ls.indexOmega,
                             indiv = TRUE,
                             name.param = pS2$name.param,
                             n.param = pS2$n.param,
                             n.cluster = pS2$n.cluster)
         
        attr(dVcov.dtheta, "score") <- out.score        
    }
    
    ## ** export
    attr(dVcov.dtheta, "vcov.param") <- vcov.param
    attr(dVcov.dtheta, "param") <- p
    return(dVcov.dtheta)      
 
}

## * dVcov2.lme
#' @rdname dVcov2
#' @export
dVcov2.lme <- dVcov2.gls

## * dVcov2.lvmfit
#' @rdname dVcov2
#' @export
dVcov2.lvmfit <- function(object, vcov.param = NULL,
                          adjust.residuals = TRUE, numericDerivative = FALSE,
                          return.score = FALSE, ...){

    detail <- NULL ## [:for CRAN check] subset
    
    p <- lava::pars(object)
    data <- stats::model.frame(object)
 
    n.param <- length(p)
    name.param <- names(p)

    as.clubSandwich <- 1

    ## ** Pre-compute quantities
    if(is.null(object$prepareScore2)){
        ## when using numerical derivative the score is computed for different sets of parameters
        ## therefore the pre-computations should not use the estimated parameters
        ## when using explicit formula for the derivative it is ok to use the estimated parameters in the pre-computations
        object$prepareScore2 <- prepareScore2(object,
                                              second.order = TRUE,
                                              usefit = (numericDerivative==FALSE) )
    }

    if(return.score || numericDerivative == FALSE){
        epsilonD2 <- residuals2(object, p = p, data = data,
                                adjust.residuals = adjust.residuals,
                                as.clubSandwich = as.clubSandwich,
                                return.prepareScore2 = TRUE, return.vcov.param = TRUE, second.order = TRUE)
        pS2  <- attr(epsilonD2, "prepareScore2")
    }
        
    ## ** param with non-zero third derivative
    keep.type <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
    tableType <- coefType(object, as.lava=FALSE)        
    keep.param <- subset(tableType,
                         subset = !is.na(lava) & detail %in% keep.type,
                         select = "originalLink", drop = TRUE)
    
    ### ** Compute the gradient 
    if(numericDerivative){

        ### *** Define function to compute the variance-covariance matrix
        if(adjust.residuals==FALSE){
            calcVcov <- function(iParam){ # x <- p.obj
                pp <- p
                pp[names(iParam)] <- iParam
                Info <- information(object, p = pp)
                return(chol2inv(chol(Info))) 
                ## return(Info) ## [:DEBUG]
            }                        
        }else{            
            calcVcov <- function(iParam){ # x <- p.obj
                pp <- p
                pp[names(iParam)] <- iParam
                res.tempo <- residuals2(object, p = pp, data = data,
                                        adjust.residuals = adjust.residuals,                                        
                                        as.clubSandwich = as.clubSandwich,
                                        return.vcov.param = TRUE)
                vcov.tempo <- attr(res.tempo, "vcov.param")
                attr(vcov.param, "warning")  <- attr(res.tempo, "warning")
                return(vcov.tempo)
            }
        }

        if(is.null(vcov.param)){
            vcov.param <- calcVcov(p)
            dimnames(vcov.param) <- list(name.param, name.param)
        }
        
        ### *** numerical derivative
        jac.param <- p[keep.param]
        res.tempo <- numDeriv::jacobian(calcVcov, x = jac.param, method = "Richardson")

        dVcov.dtheta <- array(res.tempo,
                              dim = c(n.param,n.param,length(jac.param)),
                              dimnames = list(name.param, name.param, keep.param))
        
    }else{
        vcov.param  <- attr(epsilonD2, "vcov.param")
        attr(vcov.param, "warning")  <- attr(epsilonD2, "warning")
 
        dInfo.dtheta <- .dinformation2(dmu.dtheta = pS2$dtheta$dmu.dtheta,
                                       d2mu.dtheta2 = pS2$dtheta2$d2mu.dtheta2,
                                       dOmega.dtheta = pS2$dtheta$dOmega.dtheta,
                                       d2Omega.dtheta2 = pS2$dtheta2$d2Omega.dtheta2,
                                       Omega = pS2$Omega,
                                       ls.indexOmega = pS2$ls.indexOmega,
                                       hat = pS2$hat,
                                       n.param  = pS2$n.param,
                                       name.param  = pS2$name.param,
                                       name.deriv = keep.param,
                                       n.cluster = pS2$n.cluster)

        p3 <- dim(dInfo.dtheta)[3]
        dVcov.dtheta <- array(NA, dim = dim(dInfo.dtheta), dimnames = dimnames(dInfo.dtheta))
        for(iP in 1:p3){
            dVcov.dtheta[,,iP] <- - vcov.param %*% dInfo.dtheta[,,iP] %*% vcov.param 
            ## dVcov.dtheta[,,iP] <- dInfo.dtheta[,,iP] ## [:DEBUG]
        }

    }

    ## ** score
    if(return.score){
        attr(epsilonD2, "prepareScore2") <- NULL
        attr(epsilonD2, "vcov.param") <- NULL
        out.score <- .score2(dmu.dtheta = pS2$dtheta$dmu.dtheta,
                             dOmega.dtheta = pS2$dtheta$dOmega.dtheta,
                             epsilon = epsilonD2,
                             Omega = pS2$Omega,
                             ls.indexOmega = NULL,
                             indiv = TRUE,                         
                             name.param = pS2$name.param,
                             n.param = pS2$n.param,
                             n.cluster = pS2$n.cluster)
        attr(dVcov.dtheta, "score") <- out.score
    }
  
    ## ** export
    attr(dVcov.dtheta, "vcov.param") <- vcov.param
    attr(dVcov.dtheta, "param") <- p
    return(dVcov.dtheta)       
}





## * dVcov2.lvmfit2
#' @rdname dVcov2
#' @export
dVcov2.lvmfit2 <- function(object, ...){
    class(object) <- setdiff(class(object),"lvmfit2")
    return(dVcov2(object, ...))    
}

## * dVcov2<-
#' @rdname dVcov2
#' @export
`dVcov2<-` <-
  function(x, ..., value) UseMethod("dVcov2<-")

## * dVcov2<-.lm
#' @rdname dVcov2
#' @export
`dVcov2<-.lm` <- function(x, ..., value){
    x$dVcov <- dVcov2(x, ..., adjust.residuals = value)
    class(x) <- append("lm2",class(x))

    return(x)
}    
## * dVcov2<-.gls
#' @rdname dVcov2
#' @export
`dVcov2<-.gls` <- function(x, ..., value){
    x$dVcov <- dVcov2(x, ..., adjust.residuals = value)
    class(x) <- append("gls2",class(x))

    return(x)
}    
## * dVcov2<-.lme
#' @rdname dVcov2
#' @export
`dVcov2<-.lme` <- function(x, ..., value){
    x$dVcov <- dVcov2(x, ..., adjust.residuals = value)
    class(x) <- append("lme2",class(x))
    
    return(x)
}    

## * dVcov2<-.lvmfit
#' @rdname dVcov2
#' @export
`dVcov2<-.lvmfit` <- function(x, ..., value){
    x$dVcov <- dVcov2(x, ..., adjust.residuals = value)
    class(x) <- append("lvmfit2",class(x))

    return(x)
}    

## * dVcov2<-.lvmfit2
#' @rdname dVcov2
#' @export
`dVcov2<-.lvmfit2` <- function(x, ..., value){

    class(x) <- setdiff(class(x),"lvmfit2")
    x$dVcov <- dVcov2(x, ..., adjust.residuals = value)
    class(x) <- append("lvmfit2",class(x))
    
    return(x)
}

## * .dinformation2
.dinformation2 <- function(dmu.dtheta, d2mu.dtheta2,
                           dOmega.dtheta, d2Omega.dtheta2,
                           Omega, ls.indexOmega, hat,
                           n.param, name.param, name.deriv,
                           n.cluster){

### ** prepare
    index.deriv <- match(name.deriv, name.param)
    clusterSpecific <- !is.null(ls.indexOmega)
    iOmega <- chol2inv(chol(Omega))        

    if(!clusterSpecific){
        Omega.tempo <- Omega
        iOmega.tempo <- iOmega

        ## *** small sample correction               
        df.mean <- Reduce("+",hat)
        iN.cluster <- as.double(n.cluster - diag(df.mean))
    }
    
    ### ** compute the derivative of the information matrix for each parameters
    dInfo <-  array(0,
                    dim = c(n.param, n.param, length(name.deriv)),
                    dimnames = list(name.param, name.param, name.deriv))
    
    for(iDeriv in index.deriv){ # iDeriv <- 4
        for(iP1 in 1:n.param){ # iP1 <- 1
            for(iP2 in iP1:n.param){ # iP2 <- 1
                
                iNameD <- name.param[iDeriv]
                iName1 <- name.param[iP1]
                iName2 <- name.param[iP2]

                ##cat(iNameD," ",iName1,"",iName2,"\n")
                
                test.Omega1 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dOmega.dtheta[[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2a <- !is.null(d2Omega.dtheta2[[iNameD]][[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2b <- !is.null(d2Omega.dtheta2[[iName1]][[iNameD]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega3a <- !is.null(d2Omega.dtheta2[[iNameD]][[iName2]]) && !is.null(dOmega.dtheta[[iName1]])
                test.Omega3b <- !is.null(d2Omega.dtheta2[[iName2]][[iNameD]]) && !is.null(dOmega.dtheta[[iName1]])
                
                test.mu1a <- !is.null(d2mu.dtheta2[[iNameD]][[iName1]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu1b <- !is.null(d2mu.dtheta2[[iName1]][[iNameD]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu2a <- !is.null(d2mu.dtheta2[[iNameD]][[iName2]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu2b <- !is.null(d2mu.dtheta2[[iName2]][[iNameD]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu3 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dmu.dtheta[[iName1]]) && !is.null(dmu.dtheta[[iName2]])

                ## *** Individual specific Omega (e.g. presence of missing values)
                if(clusterSpecific){
                    
                    for(iC in 1:n.cluster){

                        ## prepare
                        Omega.tempo <- Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        iOmega.tempo <- iOmega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        if(!is.null(dmu.dtheta[[iName1]])){
                            dmu.1 <- dmu.dtheta[[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dmu.dtheta[[iName2]])){
                            dmu.2 <- dmu.dtheta[[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D1 <- d2mu.dtheta2[[iNameD]][[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D1 <- d2mu.dtheta2[[iName1]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D <- d2mu.dtheta2[[iNameD]][[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D <- d2mu.dtheta2[[iName2]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iNameD]])){
                            iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName1]])){
                            iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName2]])){
                            iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega2a){
                            d2Omega.D1 <- d2Omega.dtheta2[[iNameD]][[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.Omega2b){
                            d2Omega.D1 <- d2Omega.dtheta2[[iName1]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega3a){
                            d2Omega.D2 <- d2Omega.dtheta2[[iNameD]][[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else{
                            d2Omega.D2 <- d2Omega.dtheta2[[iName2]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }

                        ## small sample correction  
                        iW.cluster <- 1 -  diag(hat[[iC]])
                        
                        ## compute
                        if(test.Omega1){                            
                            iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                            iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iW.cluster + iDiag2 * iW.cluster)
                        }
                        
                        if(test.Omega2a || test.Omega2b){
                            iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.Omega3a || test.Omega3b){
                            iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.mu1a || test.mu1b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + d2mu.D1 %*% iOmega %*% t(dmu.2)
                        }

                        if(test.mu2a || test.mu2b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + dmu.1 %*% iOmega %*% t(d2mu.D2)
                        }

                        if(test.mu3){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo %*% t(dmu.2)
                        }
                    
                }
            }
            
            ## *** Same for all individuals
                if(clusterSpecific == FALSE){

                    ## prepare
                    if(!is.null(dmu.dtheta[[iName1]])){
                        dmu.1 <- dmu.dtheta[[iName1]]
                    }
                    if(!is.null(dmu.dtheta[[iName2]])){
                        dmu.2 <- dmu.dtheta[[iName2]]
                    }
                    if(test.mu1a){
                        d2mu.D1 <- d2mu.dtheta2[[iNameD]][[iName1]]
                    }else if(test.mu1b){
                        d2mu.D1 <- d2mu.dtheta2[[iName1]][[iNameD]]
                    }
                    if(test.mu2a){
                        d2mu.D2 <- d2mu.dtheta2[[iNameD]][[iName2]]
                    }else if(test.mu2b){
                        d2mu.D2 <- d2mu.dtheta2[[iName2]][[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iNameD]])){
                        iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iName1]])){
                        iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]]
                    }
                    if(!is.null(dOmega.dtheta[[iName2]])){
                        iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]]
                    }
                    if(test.Omega2a){
                        d2Omega.D1 <- d2Omega.dtheta2[[iNameD]][[iName1]]
                    }else if(test.Omega2b){
                        d2Omega.D1 <- d2Omega.dtheta2[[iName1]][[iNameD]]
                    }
                    if(test.Omega3a){
                        d2Omega.D2 <- d2Omega.dtheta2[[iNameD]][[iName2]]
                    }else{
                        d2Omega.D2 <- d2Omega.dtheta2[[iName2]][[iNameD]]
                    }

                    ## compute
                    if(test.Omega1){
                        iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                        iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iN.cluster + iDiag2 * iN.cluster)
                    }

                    if(test.Omega2a || test.Omega2b){
                        iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.Omega3a || test.Omega3b){
                        iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.mu1a || test.mu1b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% iOmega * dmu.2)
                    }

                    if(test.mu2a || test.mu2b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu.1 %*% iOmega * d2mu.D2)
                    }

                  
                    if(test.mu3){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo * dmu.2)
                    }
            }
            
            }
        }
        dInfo[,,iNameD] <- dInfo[,,iNameD] + t(dInfo[,,iNameD]) - diag(diag(dInfo[,,iNameD]))
    }

    ### ** export
    return(dInfo)
}


##----------------------------------------------------------------------
### dVcov2.R ends here





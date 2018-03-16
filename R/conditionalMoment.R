### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: mar 16 2018 (11:46) 
##           By: Brice Ozenne
##     Update #: 930
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * conditionalMoment - documentation
#' @title Prepare the Computation of score2
#' @description Compute partial derivatives regarding to the mean and the variance, and compute the design matrices.
#' @name conditionalMoment
#' 
#' @param object,x a latent variable model.
#' @param data [data.frame] data set.
#' @param formula [formula] two-sided linear formula.
#' @param param,p [numeric vector] the fitted coefficients.
#' @param attr.param [character vector] the type of each coefficient
#' (e.g. mean or variance coefficient).
#' @param ref.group [character vector] the levels of the variable defining the variance component in a generic covariance matrix.
#' @param second.order [logical] should the terms relative to the third derivative of the likelihood be be pre-computed?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param n.cluster [integer >0] the number of i.i.d. observations.
#' @param n.endogenous [integer >0] the number of outcomes.
#' @param usefit,value [logical] If TRUE the coefficients estimated by the model are used to pre-compute quantities. Only for lvmfit objects.
#' @param name.endogenous [character vector, optional] name of the endogenous variables
#' @param name.latent [character vector, optional] name of the latent variables
#' @param ... [internal] only used by the generic method or by the <- methods.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient (\code{conditionalMoment.lvm}).
#' \item an advanced one that require the model coefficients (\code{conditionalMoment.lvmfit}). 
#' }
#' 
#' @examples
#' conditionalMoment <- lavaSearch2:::conditionalMoment
#' 
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' d <- lava::sim(m,1e2)
#' e <- estimate(m, d)
#'
#' ## basic pre-computation
#' res1 <- conditionalMoment(e, data = d, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = FALSE)
#' res1$skeleton$Sigma
#' 
#' ## full pre-computation
#' res2 <- conditionalMoment(e, param = coef(e), data = d, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = TRUE
#' )
#' res2$value$Sigma
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
`conditionalMoment` <-
  function(object, ...) UseMethod("conditionalMoment")


## * conditionalMoment.lm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lm <- function(object, name.endogenous, ...){

    X <- model.matrix(object)
    dmu <- lapply(1:NCOL(X), function(i){
        M <- X[,i,drop=FALSE]
        colnames(M) <- name.endogenous
        return(M)
    })
    names(dmu) <- colnames(X)
    
    return(list(dmu = dmu,
                dOmega = list(sigma2 = matrix(1)),
                name.3deriv = "sigma2"))
}

## * conditionalMoment.gls
#' @rdname conditionalMoment
#' @export
conditionalMoment.gls <- function(object, data, formula, 
                                  param, attr.param, ref.group,
                                  second.order,
                                  index.Omega, cluster, n.cluster, name.endogenous, 
                                  ...){
    
    ### ** prepare
    n.endogenous <- length(name.endogenous)

    ## *** coefficients
    name.varcoef <- attr.param$var.coef
    name.corcoef <- attr.param$cor.coef
    n.varcoef <- length(name.varcoef)
    n.corcoef <- length(name.corcoef)
    var.coef <- param[name.varcoef]
    cor.coef <- param[name.corcoef]

    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)

    ## *** design matrix    
    X <- stats::model.matrix(formula, data)
    X <- X[,attr.param$mean.coef,drop=FALSE] ## drop unused columns (e.g. factor with 0 occurence)    
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    
    ## *** variance terms
    if("NULL" %in% class.var == FALSE){
        name.other <- setdiff(names(var.coef),"sigma2")
        factor.varcoef <- setNames(c(1,var.coef[name.other]),
                                   attr(object$modelStruct$varStruct,"groupNames"))
        sigma2.base0 <- factor.varcoef[ref.group]        
    }else{
        sigma2.base0 <- setNames(rep(1, n.endogenous), name.endogenous)
    }
    sigma2.base <- sigma2.base0 * var.coef["sigma2"]

    ## *** corelation terms
    if("NULL" %in% class.cor == FALSE){
        M.corcoef <- matrix("", n.endogenous, n.endogenous,
                            dimnames = list(name.endogenous,name.endogenous))
        M.corcoef[which(lower.tri(M.corcoef))] <- name.corcoef
        M.corcoef <- symmetrize(M.corcoef)

        index.lower.tri <- which(lower.tri(M.corcoef))
        indexArr.lower.tri <- which(lower.tri(M.corcoef), arr.ind = TRUE)

        Msigma2.base0 <- matrix(0, n.endogenous, n.endogenous,
                                dimnames = list(name.endogenous, name.endogenous))
        Msigma2.base0[index.lower.tri] <- apply(indexArr.lower.tri, 1, function(x){sqrt(prod(sigma2.base0[x]))})
        Msigma2.base0 <- symmetrize(Msigma2.base0)
    }

### ** score - mean
    name.X <- colnames(X)
    dmu <- lapply(name.X, function(iCoef){ # iCoef <- name.X[1]
        dmu.tempo <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                            dimnames = list(NULL, name.endogenous))
        for(iC in 1:n.cluster){ ## iC <- 1
            dmu.tempo[iC,index.Omega[[iC]]] <- X[cluster==iC,iCoef]
        }
        return(dmu.tempo)
    })
    names(dmu) <- name.X

    ### ** score - variance/covariance
    dOmega <- vector(mode = "list", length = n.corcoef + n.varcoef)
    names(dOmega) <- c(name.corcoef, name.varcoef)

    ## *** dispersion coefficient
    dOmega[["sigma2"]] <- diag(sigma2.base0, nrow = n.endogenous, ncol = n.endogenous)
   
    if("NULL" %in% class.cor == FALSE){
        dOmega[["sigma2"]][index.lower.tri] <- Msigma2.base0[index.lower.tri] * cor.coef[M.corcoef[index.lower.tri]]
        dOmega[["sigma2"]] <- symmetrize(dOmega[["sigma2"]])      
    }
    dimnames(dOmega[["sigma2"]]) <-  list(name.endogenous, name.endogenous)

    ## *** other variance coefficients
    if("NULL" %in% class.var == FALSE){

        for(iVar in name.other){ # iVar <- name.other
            iTest.endogenous <- ref.group %in% iVar
            dOmega[[iVar]] <- var.coef["sigma2"]*diag(iTest.endogenous,
                                                      nrow = n.endogenous, ncol = n.endogenous)

            if("NULL" %in% class.cor == FALSE){
                index.iVar <- which(rowSums(indexArr.lower.tri==which(iTest.endogenous))>0)

                ##  d sqrt(x) / d x = 1/(2 sqrt(x)) = sqrt(x) / (2*x)
                dOmega[[iVar]][index.lower.tri[index.iVar]] <- var.coef["sigma2"]*dOmega[["sigma2"]][index.lower.tri[index.iVar]]/(2*var.coef[iVar])
                dOmega[[iVar]] <- symmetrize(dOmega[[iVar]])
            }
            
            dimnames(dOmega[[iVar]]) <- list(name.endogenous, name.endogenous)            
        }
    }
    
    ## *** correlation
    if("NULL" %in% class.cor == FALSE){
        for(iVar in name.corcoef){
            dOmega[[iVar]] <- Msigma2.base0 * var.coef["sigma2"] * (M.corcoef==iVar)
        }
    }

### ** second order
    d2Omega <- list()

    if(second.order){

        if("NULL" %in% class.var == FALSE){
            for(iVar in name.other){ ## iVar <- name.other[1]
                d2Omega[["sigma2"]][[iVar]] <- dOmega[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.cor == FALSE){
            for(iVar in name.corcoef){
                d2Omega[["sigma2"]][[iVar]] <- dOmega[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.var == FALSE && "NULL" %in% class.cor == FALSE){
            M.corvalue <- matrix(1, nrow = n.endogenous, ncol = n.endogenous)
            M.corvalue[index.lower.tri] <- cor.coef[M.corcoef[index.lower.tri]]
            M.corvalue <- symmetrize(M.corvalue, update.upper = TRUE)

            for(iVar1 in name.other){ ## iVar <- name.other[1]

                iIndex.var1 <- which(name.varcoef == iVar1)
                
                ## var var
                for(iVar2 in name.varcoef[iIndex.var1:n.varcoef]){

                    ##
                    M.tempo <- c(1,-1)[(iVar1==iVar2)+1] * dOmega[[iVar1]]/(2*var.coef[iVar2])

                    ## remove null derivative on the diagonal
                    diag(M.tempo) <- 0

                    ## remove null derivative outside the diagonal
                    iIndex.var2 <- which(name.varcoef == iVar2)
                    
                    index0 <- union(which(rowSums(indexArr.lower.tri==iIndex.var1)==0),
                                    which(rowSums(indexArr.lower.tri==iIndex.var2)==0))
                    M.tempo[index.lower.tri[index0]] <- 0
                    M.tempo <- symmetrize(M.tempo, update.upper = TRUE)

                    d2Omega[[iVar1]][[iVar2]] <- M.tempo
                }                
                
                ## var cor
                for(iVar2 in name.corcoef){                    
                    M.tempo <- dOmega[[iVar1]]/M.corvalue
                    M.tempo[M.corcoef!=iVar2] <- 0
                    if(any(M.tempo!=0)){
                        d2Omega[[iVar1]][[iVar2]] <- M.tempo
                    }
                }

            }
        }

    }

### ** export
    return(list(X = X,
                dmu = dmu,
                dOmega = dOmega,
                d2Omega = d2Omega))
    
}

## * conditionalMoment.lme
#' @rdname conditionalMoment
#' @export
conditionalMoment.lme <- function(object, attr.param, name.endogenous, ...){

    resGLS <- conditionalMoment.gls(object, attr.param = attr.param,
                                    name.endogenous = name.endogenous, ...)
        
    ### ** random effect
    name.rancoef <- attr.param$ran.coef
    n.endogenous <- length(name.endogenous)
    resGLS$dOmega[[name.rancoef]] <- matrix(1, nrow = n.endogenous, ncol = n.endogenous,
                                            dimnames = list(name.endogenous,name.endogenous)
                                            )

### ** export
    return(resGLS)
}

## * conditionalMoment.lvm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvm <- function(object, data, second.order,
                                  name.endogenous, name.latent,
                                  ...){

### ** Initialize conditional moments   
    skeleton <- skeleton(object,
                         name.endogenous = name.endogenous, 
                         name.latent = name.latent, 
                         as.lava = TRUE)
    
### ** Initialize partial derivatives of the conditional moments
    dtheta <- skeletonDtheta(object, data = data,
                             df.param.all = skeleton$df.param,
                             param2originalLink = skeleton$param2originalLink,
                             name.endogenous = name.endogenous, 
                             name.latent = name.latent)

### ** Initialize second order partial derivatives of the conditional moments
    if(second.order){
        d2theta <- skeletonDtheta2(object, data = data,
                                   df.param.all = skeleton$df.param,
                                   param2originalLink = skeleton$param2originalLink,
                                   name.latent = name.latent)
    }else{
        d2theta <- NULL
    }
    
### ** Export
    return(list(skeleton = skeleton, dtheta = dtheta, d2theta = d2theta))
}
    
    
## * conditionalMoment.lvmfit
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvmfit <- function(object, data, param, usefit, second.order, ...){

### ** normalize arguments
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)
    name.latent <- latent(object)
    n.latent <- length(name.latent)

    data <- as.matrix(data[,lava::manifest(object),drop=FALSE])

### ** initialize
    if(!is.null(object$conditionalMoment)){
        out <- object$conditionalMoment
    }else{
        out <- conditionalMoment(lava::Model(object), data = data,
                                 second.order = second.order,
                                 name.endogenous = name.endogenous, name.latent = name.latent)
    }

    ## out2 <- conditionalMoment(lava::Model(object), data = data,
    ##                              second.order = second.order,
    ##                              name.endogenous = name.endogenous, name.latent = name.latent)
### ** update according to the value of the model coefficients
    if(usefit){
        
### ** Update skeleton 
        out$skeleton <- skeleton(object, skeleton = out$skeleton,
                                 data = data, param = param,
                                 name.endogenous = name.endogenous, 
                                 name.latent = name.latent, 
                                 )
        out$skeleton$skeleton <- NULL ## cleaning
 
        ### ** Update first order partial derivatives
        out$dtheta <- skeletonDtheta(object, dtheta = out$dtheta,
                                     name.endogenous = name.endogenous, 
                                     name.latent = name.latent,
                                     B = out$skeleton$value$B,
                                     alpha.XGamma = out$skeleton$value$alpha.XGamma,
                                     Lambda = out$skeleton$value$Lambda,
                                     Psi = out$skeleton$value$Psi)

        ### ** Compute second order partial derivatives
        if(second.order){
            out$d2theta <- skeletonDtheta2(object,
                                           dtheta = out$dtheta,
                                           d2theta = out$d2theta,
                                           name.endogenous = name.endogenous,
                                           name.latent = name.latent,
                                           B = out$skeleton$value$B,
                                           Lambda = out$skeleton$value$Lambda,
                                           Psi = out$skeleton$value$Psi)
            out$d2theta <- out$d2theta[c("d2mu","d2Omega")]
        }
        out$skeleton$value$alpha.XGamma.iIB <- out$dtheta$alpha.XGamma.iIB
        out$skeleton$value$tLambda.tiIB.Psi.iIB <- out$dtheta$tLambda.tiIB.Psi.iIB
        out$dtheta <- out$dtheta[c("dmu","dOmega")]
    }
    
    ### ** Export
    return(c(out[["skeleton"]], out[["dtheta"]], out[["d2theta"]]))    
}


### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: feb  5 2018 (17:19) 
##           By: Brice Ozenne
##     Update #: 766
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * prepareScore2 - documentation
#' @title Prepare the Computation of score2
#' @description Compute partial derivatives regarding to the mean and the variance, and compute the design matrices.
#' @name prepareScore2
#' 
#' @param object,x a latent variable model.
#' @param X [matrix] the design matrix.
#' @param param,p [numeric vector] the fitted coefficients.
#' @param attr.param [character vector] the type of each coefficient
#' (e.g. mean or variance coefficient).
#' @param second.order [logical] should the terms relative to the third derivative of the likelihood be be pre-computed?
#' @param n.cluster [integer >0] the number of i.i.d. observations.
#' @param n.endogenous [integer >0] the number of outcomes.
#' @param index.obs [integer vector] the indexes of the outcomes relative to each observation (e.g. 1,3 if only outcome 1 and 3 are observed for the observation).
#' @param usefit,value [logical] If TRUE the coefficients estimated by the model are used to pre-compute quantities. Only for lvmfit objects.
#' @param data [data.frame, optional] data set.
#' @param name.endogenous [character vector, optional] name of the endogenous variables
#' @param name.latent [character vector, optional] name of the latent variables
#' @param ... [internal] only used by the generic method or by the <- methods.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient
#' \item an advanced one that require the model coefficients. 
#' }
#' 
#' @examples
#' prepareScore2 <- lavaSearch2:::prepareScore2
#' 
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m, lava::sim(m,1e2))
#' res <- prepareScore2(e)
#' res$skeleton$df.param
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
`prepareScore2` <-
  function(object, ...) UseMethod("prepareScore2")

## * prepareScore2.gls
#' @rdname prepareScore2
#' @export
prepareScore2.gls <- function(object, X, 
                              param, attr.param,
                              second.order,
                              n.cluster, n.endogenous, name.endogenous, index.obs,
                              ...){
    
### ** prepare
    
    ## *** coefficients
    name.varcoef <- attr.param$var.coef
    name.corcoef <- attr.param$cor.coef
    n.varcoef <- length(name.varcoef)
    n.corcoef <- length(name.corcoef)
    var.coef <- param[name.varcoef]
    cor.coef <- param[name.corcoef]

    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)

    ## *** variance terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if("NULL" %in% class.var){            
        sigma2.base <- stats::setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }else{
        sigma2.base <- stats::setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }
    sigma2.base0 <- sigma2.base/var.coef["sigma2"]

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
    dmu.dtheta <- lapply(name.X, function(iCoef){
        M.tempo <- matrix(NA, nrow = n.cluster, ncol = n.endogenous)    
        M.tempo[index.obs] <- X[,iCoef]
        colnames(M.tempo) <- name.endogenous
        return(M.tempo)
    })
    names(dmu.dtheta) <- name.X

### ** score - variance/covariance
    dOmega.dtheta <- vector(mode = "list", length = n.corcoef + n.varcoef)
    names(dOmega.dtheta) <- c(name.corcoef, name.varcoef)
    
    ## *** dispersion coefficient
    dOmega.dtheta[["sigma2"]] <- diag(sigma2.base0[name.endogenous], nrow = n.endogenous, ncol = n.endogenous)
   
    if("NULL" %in% class.cor == FALSE){
        dOmega.dtheta[["sigma2"]][index.lower.tri] <- Msigma2.base0[index.lower.tri] * cor.coef[M.corcoef[index.lower.tri]]
        dOmega.dtheta[["sigma2"]] <- symmetrize(dOmega.dtheta[["sigma2"]])      
    }
    dimnames(dOmega.dtheta[["sigma2"]]) <-  list(name.endogenous, name.endogenous)

    ## *** other variance coefficients
    if("NULL" %in% class.var == FALSE){
        for(iVar in name.other){
            dOmega.dtheta[[iVar]] <- var.coef["sigma2"]*diag(name.endogenous %in% iVar,
                                                             nrow = n.endogenous, ncol = n.endogenous)

            if("NULL" %in% class.cor == FALSE){
                iEndogenous <- which(name.endogenous==iVar)
                index.iVar <- which(rowSums(indexArr.lower.tri==iEndogenous)>0)

                ##  d sqrt(x) / d x = 1/(2 sqrt(x)) = sqrt(x) / (2*x)
                dOmega.dtheta[[iVar]][index.lower.tri[index.iVar]] <- var.coef["sigma2"]*dOmega.dtheta[["sigma2"]][index.lower.tri[index.iVar]]/(2*var.coef[iVar])
                dOmega.dtheta[[iVar]] <- symmetrize(dOmega.dtheta[[iVar]])
            }
            
            dimnames(dOmega.dtheta[[iVar]]) <- list(name.endogenous, name.endogenous)            
        }
    }
    
    ## *** correlation
    if("NULL" %in% class.cor == FALSE){
        for(iVar in name.corcoef){
            dOmega.dtheta[[iVar]] <- Msigma2.base0 * var.coef["sigma2"] * (M.corcoef==iVar)
        }
    }

### ** second order
    d2Omega.dtheta2 <- list()

    if(second.order){

        if("NULL" %in% class.var == FALSE){
            for(iVar in name.other){ ## iVar <- name.other[1]
                d2Omega.dtheta2[["sigma2"]][[iVar]] <- dOmega.dtheta[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.cor == FALSE){
            for(iVar in name.corcoef){
                d2Omega.dtheta2[["sigma2"]][[iVar]] <- dOmega.dtheta[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.var == FALSE && "NULL" %in% class.cor == FALSE){
            M.corvalue <- matrix(1, nrow = n.endogenous, ncol = n.endogenous)
            M.corvalue[index.lower.tri] <- cor.coef[M.corcoef[index.lower.tri]]
            M.corvalue <- symmetrize(M.corvalue, update.upper = TRUE)

            for(iVar1 in name.other){ ## iVar <- name.other[1]

                iIndex.var1 <- which(names(sigma2.base0) == iVar1)
                
                ## var var
                for(iVar2 in name.varcoef[iIndex.var1:n.varcoef]){

                    ##
                    M.tempo <- c(1,-1)[(iVar1==iVar2)+1] * dOmega.dtheta[[iVar1]]/(2*var.coef[iVar2])

                    ## remove null derivative on the diagonal
                    diag(M.tempo) <- 0

                    ## remove null derivative outside the diagonal
                    iIndex.var2 <- which(name.varcoef == iVar2)
                    
                    index0 <- union(which(rowSums(indexArr.lower.tri==iIndex.var1)==0),
                                    which(rowSums(indexArr.lower.tri==iIndex.var2)==0))
                    M.tempo[index.lower.tri[index0]] <- 0
                    M.tempo <- symmetrize(M.tempo, update.upper = TRUE)

                    ##
                    ##if(iVar1==iVar2){
                    d2Omega.dtheta2[[iVar1]][[iVar2]] <- M.tempo
                    ##}
                }                
                
                ## var cor
                for(iVar2 in name.corcoef){                    
                    M.tempo <- dOmega.dtheta[[iVar1]]/M.corvalue
                    M.tempo[M.corcoef!=iVar2] <- 0
                    if(any(M.tempo!=0)){
                        d2Omega.dtheta2[[iVar1]][[iVar2]] <- M.tempo
                    }
                }

            }
        }

    }
    
### ** export
    return(list(dmu.dtheta = dmu.dtheta,
                dOmega.dtheta = dOmega.dtheta,
                d2Omega.dtheta2 = d2Omega.dtheta2))
    
}

## * prepareScore2.lme
#' @rdname prepareScore2
#' @export
prepareScore2.lme <- function(object, X, 
                              param, attr.param,
                              second.order,
                              n.cluster, n.endogenous, name.endogenous, index.obs, ...){

    resGLS <- prepareScore2.gls(object, X = X, param = param, attr.param = attr.param,
                                n.cluster = n.cluster,
                                second.order = second.order,
                                n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                index.obs = index.obs)
        
### ** random effect
    name.rancoef <- attr.param$ran.coef
    resGLS$dOmega.dtheta[[name.rancoef]] <- matrix(1, nrow = n.endogenous, ncol = n.endogenous,
                                                   dimnames = list(name.endogenous,name.endogenous))

### ** export
    return(resGLS)
}
    

#' @rdname prepareScore2
#' @export
`prepareScore2<-` <-
  function(object, ..., value) UseMethod("prepareScore2<-")

#' @rdname prepareScore2
#' @export
"prepareScore2<-.lvmfit" <- function(object, ..., value) {
    object$prepareScore2  <- prepareScore2(lava::Model(object), data = value, ...)
    return(object)
}
## * prepareScore2.lvm
#' @rdname prepareScore2
#' @export
prepareScore2.lvm <- function(object, data, second.order,
                              name.endogenous = NULL, name.latent = NULL,
                              ...){

    pS2 <- list()
    
    ### ** Compute skeleton   
    pS2$skeleton <- skeleton(object,
                                       name.endogenous = name.endogenous, 
                                       name.latent = name.latent, 
                                       as.lava = TRUE)
    
    ### ** Initialize partial derivatives
    pS2$dtheta <- skeletonDtheta(object, data = data,
                                           df.param.all = pS2$skeleton$df.param,
                                           param2originalLink = pS2$skeleton$param2originalLink,
                                           name.endogenous = name.endogenous, 
                                           name.latent = name.latent)

    ### ** Initialize second order partial derivatives
    if(second.order){
        pS2$dtheta2 <- skeletonDtheta2(object, data = data,
                                                 df.param.all = pS2$skeleton$df.param,
                                                 param2originalLink = pS2$skeleton$param2originalLink,
                                                 name.latent = name.latent)
    }
    
    ### ** Export
    pS2$toUpdate <- TRUE
    return(pS2)
}
    
    
## * prepareScore2.lvmfit
#' @rdname prepareScore2
#' @export
prepareScore2.lvmfit <- function(object, data = NULL, p = NULL, usefit = TRUE,
                                 name.endogenous = NULL, name.latent = NULL,
                                 second.order = FALSE, ...){

    ### ** normalize arguments
    if(is.null(name.endogenous)){name.endogenous <- endogenous(object)}
    n.endogenous <- length(name.endogenous)
    if(is.null(name.latent)){name.latent <- latent(object)}
    n.latent <- length(name.latent)

    if(is.null(data)){
        data <- stats::model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
    if(usefit==FALSE){
        pS2 <- prepareScore2(lava::Model(object), data = data, second.order = second.order,
                             name.endogenous = name.endogenous, name.latent = name.latent)
        return(pS2)    
    }
    
    if(is.null(p)){
        p <- pars(object)        
    }
    
    pS2 <- list()
    ### ** Update skeleton with current estimates
    pS2$skeleton <- skeleton(object, data = data, p = p,
                             name.endogenous = name.endogenous, 
                             name.latent = name.latent, 
                             as.lava = TRUE)
 
    ### ** Update first order partial derivatives with current estimates
    pS2$dtheta <- skeletonDtheta(object, data = data,
                                 df.param.all = pS2$skeleton$df.param,
                                 param2originalLink = pS2$skeleton$param2originalLink,
                                 name.endogenous = name.endogenous, 
                                 name.latent = name.latent,
                                 B = pS2$skeleton$value$B,
                                 alpha.XGamma = pS2$skeleton$value$alpha.XGamma,
                                 Lambda = pS2$skeleton$value$Lambda,
                                 Psi = pS2$skeleton$value$Psi)

    ### ** Compute second order partial derivatives with current estimates
    if(second.order){
        pS2$dtheta2 <- skeletonDtheta2(object, data = data,
                                       OD = pS2$dtheta,
                                       df.param.all = pS2$skeleton$df.param,
                                       param2originalLink = pS2$skeleton$param2originalLink,
                                       name.endogenous = name.endogenous,
                                       name.latent = name.latent,
                                       B = pS2$skeleton$value$B,
                                       Lambda = pS2$skeleton$value$Lambda,
                                       Psi = pS2$skeleton$value$Psi)
    }
    
    ### ** Export
    pS2$toUpdate <- FALSE
    return(pS2)    
    
}

## * prepareScore2<-
#' @rdname prepareScore2
#' @export
`prepareScore2<-` <-
  function(x, ..., value) UseMethod("prepareScore2<-")

## * prepareScore2<-.lvmfit
#' @rdname prepareScore2
#' @export
`prepareScore2<-.lvmfit` <- function(x, ..., value){
    x$prepareScore2 <- prepareScore2(x, ..., usefit = value)
    return(x)
}    




#----------------------------------------------------------------------
### prepareScore2.R ends here


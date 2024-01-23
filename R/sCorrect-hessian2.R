### sCorrect-hessian2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: jan 23 2024 (10:25) 
##           By: Brice Ozenne
##     Update #: 144
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - hessian2
#' @title  Hessian With Small Sample Correction.
#' @description  Extract the hessian from a latent variable model, with small sample correction
#' @name hessian2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param indiv [logical] If \code{TRUE}, the hessian relative to each observation is returned. Otherwise the total hessian is returned.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the hessian.
#'
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#' 
#' @return An array containing the second derivative of the likelihood relative to each sample (dim 3)
#' and each pair of model coefficients (dim 1,2).
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
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' hessian2(e.lvm)
#'
#' @concept small sample inference
#' @export
`hessian2` <-
  function(object, indiv, cluster, as.lava, ...) UseMethod("hessian2")

## * hessian2.lvmfit
#' @rdname hessian2
#' @export
hessian2.lvmfit <- function(object, indiv = FALSE, cluster = NULL, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(hessian2(estimate2(object, ssc = ssc, hessian = TRUE, ...), cluster = cluster, as.lava = as.lava))

}


## * hessian2.lvmfit2
#' @rdname hessian2
#' @export
hessian2.lvmfit2 <- function(object, indiv = FALSE, cluster = NULL, as.lava = TRUE, ...){
    
    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    ## ** define cluster
    if(is.null(cluster)){
        n.cluster <- object$sCorrect$cluster$n.cluster
        cluster.index <- 1:n.cluster
    }else{
        if(!is.numeric(cluster)){
            data <- object$sCorrect$data
            if(length(cluster)==1){                
                if(cluster %in% names(data) == FALSE){
                    stop("Invalid \'cluster\' argument \n",
                         "Could not find variable \"",cluster,"\" in argument \'data\' \n")
                }
                cluster <- data[[cluster]]
            }
            cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))            
        }else{
            cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))
        }

        n.cluster <- length(unique(cluster.index))
    }

    ## ** get hessian
    hessian <- object$sCorrect$hessian
    if(is.null(hessian)){
        hessian <- .hessian2(dmu = object$sCorrect$dmoment$dmu,
                             dOmega = object$sCorrect$dmoment$dOmega,
                             d2mu = object$sCorrect$d2moment$d2mu,
                             d2Omega = object$sCorrect$d2moment$d2Omega,
                             epsilon = object$sCorrect$residuals,                                     
                             OmegaM1 = object$sCorrect$moment$OmegaM1.missing.pattern,
                             missing.pattern = object$sCorrect$missing$pattern,
                             unique.pattern = object$sCorrect$missing$unique.pattern,
                             name.pattern = object$sCorrect$missing$name.pattern,
                             grid.mean = object$sCorrect$skeleton$grid.dmoment$mean, 
                             grid.var = object$sCorrect$skeleton$grid.dmoment$var, 
                             grid.hybrid = object$sCorrect$skeleton$grid.dmoment$hybrid, 
                             name.param = object$sCorrect$skeleton$Uparam,
                             leverage = object$sCorrect$leverage,
                             n.cluster = object$sCorrect$cluster$n.cluster,
                             weights = object$sCorrect$weights)

        hessian.name <- stats::setNames(names(object$sCorrect$skeleton$originalLink2param),object$sCorrect$skeleton$originalLink2param)[object$sCorrect$skeleton$Uparam]
        dimnames(hessian) <- list(as.character(hessian.name),
                                  as.character(hessian.name),
                                  NULL)
    }

    if(!is.null(cluster)){ ## aggregate hessian by cluster
        hessianSave <- hessian
        hessian <- array(0, dim = dim(hessian),
                         dimnames = dimnames(hessian))
        for(i in 1:length(cluster)){
            hessian[,,cluster[i]] <- hessian[,,cluster[i]] + hessianSave[,,i]
        }
    }

    ## ** export
    hessian <- hessian[names(object$sCorrect$skeleton$originalLink2param),
                       names(object$sCorrect$skeleton$originalLink2param),
                      ,
                       drop=FALSE]
    if(as.lava == FALSE){
        dimnames(hessian) <- list(as.character(object$sCorrect$skeleton$originalLink2param),
                                  as.character(object$sCorrect$skeleton$originalLink2param),
                                  NULL)
    }
    if(indiv==FALSE){
        hessian <- apply(hessian, 1:2, sum)
        
    }else if(!is.null(cluster)){
        index2.cluster <- tapply(1:length(cluster),cluster,list)
        attr(hessian,"cluster") <- names(index2.cluster)
    }
    
    return(hessian)
}

## * .hessian2
#' @title Compute the Hessian Matrix From the Conditional Moments
#' @description Compute the Hessian matrix from the conditional moments.
#' @name hessian2-internal
#' 
#' @details \code{calc_hessian} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.hessian2 <- function(dmu, dOmega, d2mu, d2Omega, epsilon, OmegaM1,
                      missing.pattern, unique.pattern, name.pattern,
                      grid.mean, grid.var, grid.hybrid, name.param,
                      leverage, n.cluster, weights){
    if(lava.options()$debug){cat(".hessian2\n")}

    ## ** Prepare
    n.grid.mean <- NROW(grid.mean)
    n.grid.var <- NROW(grid.var)
    n.grid.hybrid <- NROW(grid.hybrid)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    hessian <- array(NA, dim = c(n.param, n.param, n.cluster),
                     dimnames = list(name.param,name.param,NULL))
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
    if(length(dmu)>0 && length(dOmega)>0){
        index.hybrid <- 1:n.grid.hybrid
    }else{
        index.hybrid <- NULL
    } 

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        iEpsilon <- epsilon[iIndex,iY,drop=FALSE]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)
        if(!is.null(leverage)){
            iLeverage <- leverage[iIndex,iY,drop=FALSE]
        }
        hessian[,,iIndex] <- 0 ##  ## initialize (keep NA for missing values)

        ## *** second derivative relative to the mean parameters
        for(iG in index.mean){ # iG <- 1
            iP1 <- grid.mean[iG,1]
            iP2 <- grid.mean[iG,2]
            if(grid.mean[iG,"d2.12"]){
                term1 <- rowSums((id2mu[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.mean[iG,"d2.21"]){
                term1 <- rowSums((id2mu[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            term2 <- -rowSums((idmu[[iP1]] %*% iOmegaM1) * idmu[[iP2]])
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        ## *** second derivative relative to the variance parameters
        for(iG in index.var){ # iG <- 2
            iP1 <- grid.var[iG,1]
            iP2 <- grid.var[iG,2]

            term1a <- - diag(iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]])
            term2 <- - rowSums((iEpsilon %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            if(grid.var[iG,"d2.12"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP1]][[iP2]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.var[iG,"d2.21"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP2]][[iP1]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1b <- 0
                term3 <- 0
            }
            if(is.null(leverage)){
                hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] - 1/2 * rep(sum(term1a + term1b), length(iIndex)) + term2 + term3
            }else{
                hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] - 1/2 * rowSums( sweep(1-iLeverage, FUN = "*", STATS = term1a + term1b, MARGIN = 2) ) + term2 + term3
            }
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        
        ## *** second derivative relative to the mean and variance parameters
        for(iG in index.hybrid){ # iG <- 1
            iP1 <- grid.hybrid[iG,1]
            iP2 <- grid.hybrid[iG,2]

            if(!is.null(idmu[[iP1]]) && !is.null(idOmega[[iP2]])){
                term1 <- - rowSums((idmu[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            if(!is.null(idmu[[iP2]]) && !is.null(idOmega[[iP1]])){
                term2 <- - rowSums((idmu[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term2 <- 0
            }
            
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
    }

    ## ** export
    if(!is.null(weights)){
        for(iI in 1:length(weights)){
            hessian[,,iI] <- weights[iI] * hessian[,,iI]
        }
    }
    return(hessian)
}

## * .subsetList
.subsetList <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = function(iL){iL[indexRow,indexCol,drop=FALSE]}))
    }
}
## * .subsetList2
.subsetList2 <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = .subsetList, indexRow = indexRow, indexCol = indexCol))
    }
}


######################################################################
### sCorrect-hessian2.R ends here

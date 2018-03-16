### adjustEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: mar 15 2018 (17:37) 
##           By: Brice Ozenne
##     Update #: 502
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * adjustEstimate
#' @title Compute Bias Corrected Quantities.
#' @description Compute bias corrected residuals variance covariance matrix
#' and information matrix.
#' Also provides the leverage values and corrected sample size when adjust.n is set to TRUE.
#' 
#' @keywords internal
adjustEstimate <- function(epsilon, Omega, dmu, dOmega, n.cluster,
                           name.param, name.endogenous, name.meanparam, name.varparam,
                           index.Omega,
                           adjust.Omega, adjust.n, tol, n.iter, trace){

    ## ** Prepare
    name.hybridparam <- intersect(name.meanparam, name.varparam)

    n.param <- length(name.param)
    n.meanparam <- length(name.meanparam)
    n.varparam <- length(name.varparam)
    n.hybridparam <- length(name.hybridparam)

    n.endogenous <- length(name.endogenous)
    grid.meanparam <- .combination(name.meanparam, name.meanparam)    
    n.grid.meanparam <- NROW(grid.meanparam)
    grid.varparam <- .combination(name.varparam, name.varparam)
    n.grid.varparam <- NROW(grid.varparam)

    ## check low diagonal
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.meanparam[,1]]-name2num[grid.meanparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (mean parameter) \n")
    }
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.varparam[,1]]-name2num[grid.varparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (variance parameter) \n")
    }
    ##
    
    param2index <- setNames(1:n.param, name.param)

    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))
    ls.dmu <- vector(mode = "list", length = n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        if(is.null(index.Omega)){            
            leverage[iC,] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = n.endogenous,
                                   dimnames = list(name.param, name.endogenous))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[iC,]}))
        }else{
            leverage[iC,index.Omega[[iC]]] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = length(index.Omega[[iC]]),
                                   dimnames = list(name.param, name.endogenous[index.Omega[[iC]]]))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[iC,index.Omega[[iC]]]}))
        }        
    }
    ## ** Initialisation (i.e. first iteration without correction)
    if(any(eigen(Omega)$value<=0)){
        stop("the residual variance-covariance matrix is not positive definite \n")
    }

    if(is.null(index.Omega)){
        n.corrected <- rep(n.cluster, n.endogenous)
    }else{
        n.corrected <- NULL
    }
    ls.Psi <- vector(mode = "list", length = n.cluster)

    Omega.adj <- Omega
    if(adjust.n){
        epsilon.adj <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                              dimnames = list(NULL, name.endogenous))
    }else{
        epsilon.adj <- epsilon
    }

    if(trace>0){
        cat("* Reconstruct estimated information matrix ")
    }

    iInfo <- .information2(dmu = dmu,
                           dOmega = dOmega,
                           Omega = Omega,
                           n.corrected = n.corrected,
                           leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                           grid.meanparam = grid.meanparam,
                           n.grid.meanparam = n.grid.meanparam,
                           grid.varparam = grid.varparam,
                           n.grid.varparam = n.grid.varparam,
                           name.param = name.param,
                           name.meanparam = name.meanparam,
                           name.varparam = name.varparam,
                           param2index = param2index, n.param = n.param)
    iVcov.param <- chol2inv(chol(iInfo))
    if(trace>0){
        cat("- done \n")
    }
    
    ## ** Loop    
    if(adjust.Omega || adjust.n){
        if(trace>0){
            cat("* iterative small sample correction: ")
        }
        iIter <- 0
        iTol <- Inf
        Omega_save <- Omega     
    }else{
        iIter <- Inf
        iTol <- -Inf        
    }
    
    while(iIter < n.iter & iTol > tol){
        if(trace>0){
            cat("*")
        }

        ## *** Step (i-ii): compute individual bias, expected bias
        Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                      dimnames = list(name.endogenous, name.endogenous))
        M.countCluster <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                                 dimnames = list(name.endogenous, name.endogenous))
        for(iC in 1:n.cluster){
            ## individual bias
            ls.Psi[[iC]] <- t(ls.dmu[[iC]])  %*% iVcov.param %*% ls.dmu[[iC]]
            ## cumulated bias            
            if(is.null(index.Omega)){
                Psi <- Psi + ls.Psi[[iC]]
                M.countCluster <- M.countCluster + 1
            }else{
                Psi[index.Omega[[iC]],index.Omega[[iC]]] <- Psi[index.Omega[[iC]],index.Omega[[iC]]] + ls.Psi[[iC]]
                M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] <- M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] + 1
            }
        }

        ## update
        for(iPsi in 1:length(Psi)){
            if(M.countCluster[iPsi]>0){
                Psi[iPsi] <- Psi[iPsi]/M.countCluster[iPsi]
            }
        }
        
        ## *** Step (iii): corrected sample size
        if(adjust.n){
            ## symmetric square root.
            Omega.adj.chol <- matrixPower(Omega.adj, symmetric = TRUE, power = 1/2)
            ## Omega.adj.chol %*% Omega.adj.chol - Omega.adj
            ## Omega.adj.chol %*% Omega.adj %*% Omega.adj.chol - Omega.adj %*% Omega.adj
            if(is.null(index.Omega)){
                iOmega.adj <- Omega.adj
                iOmega.adj.chol <- matrixPower(iOmega.adj, symmetric = TRUE, power = 1/2)
                iOmegaM1.adj <- chol2inv(iOmega.adj.chol)
                iIndex.Omega <- 1:n.endogenous
            }
            
            for(iC in 1:n.cluster){                 # iC <- 1
                if(!is.null(index.Omega)){
                    iIndex.Omega <- index.Omega[[iC]]
                    iOmega.adj <- Omega.adj[iIndex.Omega,iIndex.Omega,drop=FALSE]
                    iOmega.adj.chol <- matrixPower(iOmega.adj, symmetric = TRUE, power = 1/2)
                    iOmegaM1.adj <- chol2inv(iOmega.adj.chol)
                }
                ## corrected epsilon
                H <- iOmega.adj %*% iOmega.adj - iOmega.adj.chol %*% ls.Psi[[iC]] %*% iOmega.adj.chol
                ## iH <- matrixPower(H, symmetric = TRUE, power = -1/2)
                iH <- tryCatch(matrixPower(H, symmetric = TRUE, power = -1/2), warning = function(w){w})
                if(inherits(iH,"warning")){
                    stop("Cannot compute the adjusted residuals \n",
                         "Estimated bias too large compared to the estimated variance-covariance matrix \n",
                         "Consider setting argument \'adjust.n\' to FALSE when calling sCorrect \n")
                }
                ## Omega.adj.chol %*% iH %*% Omega.adj.chol %*% (iOmega.adj - ls.Psi[[iC]]) %*% Omega.adj.chol %*% iH %*% Omega.adj.chol - Omega.adj
                epsilon.adj[iC,iIndex.Omega] <- epsilon[iC,iIndex.Omega] %*% iOmega.adj.chol %*% iH %*% iOmega.adj.chol
                ## derivative of the score regarding Y
                scoreY <- ls.dmu[[iC]] %*% iOmegaM1.adj
                for(iP in 1:n.varparam){ ## iP <- 1
                    scoreY[name.varparam[iP],] <- scoreY[name.varparam[iP],] + 2 * epsilon.adj[iC,iIndex.Omega] %*% (iOmegaM1.adj %*% dOmega[[name.varparam[iP]]][iIndex.Omega,iIndex.Omega]  %*% iOmegaM1.adj)
                }
                ## leverage
                leverage[iC,iIndex.Omega] <- colSums(iVcov.param %*% ls.dmu[[iC]] * scoreY) ## NOTE: dimensions of ls.dmu and scoreY matches even when there are missing values
                                        # same as
                                        # diag(t(ls.dmu[[iC]])  %*% iVcov.param %*% scoreY)
            }
            ## corrected sample size                
            n.corrected <- rep(n.cluster, n.endogenous) - colSums(leverage, na.rm = TRUE)
        }

        ## *** Step (iv): correct residual covariance matrix
        if(adjust.Omega){
            ## corrected residual covariance variance
            Omega.adj <- Omega + Psi
        }

        ## *** Step (v): expected information matrix
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega.adj,
                               n.corrected = n.corrected,
                               leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                               grid.meanparam = grid.meanparam,
                               n.grid.meanparam = n.grid.meanparam,
                               grid.varparam = grid.varparam,
                               n.grid.varparam = n.grid.varparam,
                               name.param = name.param,
                               name.meanparam = name.meanparam,
                               name.varparam = name.varparam,
                               param2index = param2index, n.param = n.param)
        iVcov.param <- chol2inv(chol(iInfo))
        
        ## *** Update cv
        iIter <- iIter + 1
        iTol <- norm(Omega.adj-Omega_save, type = "F")
        Omega_save <- Omega.adj
        ## cat("Omega.adj: ",Omega.adj," | n:",n.corrected," | iTol:",iTol,"\n")
    }
    
    ## ** Post processing
    if(!is.infinite(iIter)){

        if(iTol > tol){
            warning("small sample correction did not reach convergence after ",iIter," iterations \n")

            if(trace>0){
                cat(" - incomplete \n")
            }
        }else{
            if(trace>0){
                cat(" - done \n")
            }
        }
        
    }

    vcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if("try-error" %in% class(vcov.param)){
        errorMessage <- vcov.param
        vcov.param <- solve(iInfo)
        attr(vcov.param, "warning") <- errorMessage
    }
    dimnames(vcov.param) <- dimnames(iInfo)

    if(!is.null(index.Omega)){
        n.corrected <- colSums(leverage, na.rm = TRUE)
        OmegaM1 <- lapply(1:n.cluster, function(iC){
            return(solve(Omega.adj[index.Omega[[iC]],index.Omega[[iC]]]))
        })    
    }else{
        OmegaM1 <- chol2inv(chol(Omega.adj))
    }

    ## ** Export    
    return(list(Omega = Omega.adj,
                OmegaM1 = OmegaM1,
                vcov.param = vcov.param,
                leverage = leverage,
                n.corrected = n.corrected,
                epsilon = epsilon.adj,
                iter = iIter))
}


##----------------------------------------------------------------------
### adjustEstimate.R ends here

### sCorrect-sscResiduals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: jan 17 2022 (14:06) 
##           By: Brice Ozenne
##     Update #: 1226
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .init_sscResiduals
.init_sscResiduals <- function(object){

    out <- list()

    ## ** extract info
    endogenous <- object$sCorrect$endogenous
    n.endogenous <- length(endogenous)
    latent <- object$sCorrect$latent
    n.latent <- length(latent)
    
    type <- object$sCorrect$skeleton$type

    index.var <- which(type$detail %in% c("Sigma_var","Sigma_cov","Psi_var","Psi_cov"))
    index.param <- which(!is.na(type$param))
    type.param <- type[intersect(index.var,index.param),,drop=FALSE]
    
    
    Omega <- object$sCorrect$moment$Omega
    ## name.var <- object$sCorrect$name.param[object$sCorrect$name.param %in% unique(type.param$param)] ## make sure to keep the same order as in the original vector of parameters
    name.var <- unique(type.param$param)
    n.var <- length(name.var)

    ## ** subset residual variance-covariance
    index.upper.tri <- data.frame(index = which(upper.tri(Omega, diag = TRUE)),
                                  which(upper.tri(Omega, diag = TRUE), arr.ind = TRUE)
                                  )
    name.rhs <- paste(endogenous[index.upper.tri[,"row"]],
                      lava.options()$symbols[2],
                      endogenous[index.upper.tri[,"col"]],
                      sep = "")
    n.rhs <- length(name.rhs)

    ## ** fixed part of the variance-covariance matrix
    index.value <- which(!is.na(type$value))
    index.var2 <- which(type$detail %in% c("Lambda","B","Sigma_var","Sigma_cov","Psi_var","Psi_cov"))
    type.var.constrain <- type[intersect(index.value, index.var2),,drop=FALSE]

    Omega.constrain <- matrix(0, nrow = n.endogenous, ncol = n.endogenous, 
                              dimnames = list(endogenous,endogenous))

    if(NROW(type.var.constrain)>0){
        ## extract parameters with value fixed by the user
        value <- object$sCorrect$skeleton$value

        ## add fixed Sigma value
        if(any(c("Sigma_var","Sigma_cov") %in% type.var.constrain$detail)){
            addSigma <- value$Sigma
            addSigma[is.na(addSigma)] <- 0
            Omega.constrain <- Omega.constrain + addSigma
        }
        if(any(c("Psi_var","Psi_cov") %in% type.var.constrain$detail)){
            Psi.constrain <- value$Psi
            Lambda.constrain <- value$Lambda
            if(any(is.na(value$B))){
                stop("Current implementation cannot handle constraints in Psi when there are B parameters.\n")
            }
            if(any(is.na(value$Lambda))){
                stop("Current implementation cannot handle constraints in Psi when there are lambda parameters.\n")
            }
            if("B" %in% names(value)){
                iIB.constrain <- solve(diag(1, nrow = n.latent, ncol = n.latent) - value$B)
            }else{
                iIB.constrain <- diag(1, nrow = n.latent, ncol = n.latent)
            }

            addPsi <- t(Lambda.constrain) %*% t(iIB.constrain) %*% Psi.constrain %*% iIB.constrain %*% Lambda.constrain
            addPsi[is.na(addPsi)] <- 0
            Omega.constrain <- Omega.constrain + addPsi
        }

    }

    ## ** Design matrix
    A <- matrix(0, nrow = n.rhs, ncol = n.var,
                dimnames = list(name.rhs, name.var))

    ## *** Sigma_var and Sigma_cov
    if(any(type.param$detail %in% c("Sigma_var","Sigma_cov"))){
        type.Sigma <- type.param[type.param$detail %in% c("Sigma_var","Sigma_cov"),,drop=FALSE]
        for(iRow in 1:NROW(type.Sigma)){ ## iRow <- 1 
            A[paste0(type.Sigma[iRow,"Y"],"~~",type.Sigma[iRow,"X"]),type.Sigma[iRow,"param"]] <- 1
        }
    }
    attr(A,"name") <- name.var

    ## *** Psi_var and Psi_cov
    if(any(type.param$detail %in% c("Psi_var","Psi_cov"))){
        type.Psi <- type.param[type.param$detail %in% c("Psi_var","Psi_cov"),,drop=FALSE]
        index.Psi <- cbind(row = match(type.Psi$X, latent),
                           col = match(type.Psi$Y, latent))
        rownames(index.Psi) <- type.Psi$param
    }else{
        index.Psi <- NULL
    }
    
    return(list(type = "residuals",
                param0 = object$sCorrect$param,
                Omega0 = object$sCorrect$moment$Omega,
                residuals0 = object$sCorrect$residuals,
                index.upper.tri = index.upper.tri,
                name.rhs = name.rhs,
                name.var = name.var,
                A = A,
                Omega.constrain = Omega.constrain,
                index.Psi = index.Psi
                ))
}

## * .sscResiduals
#' @title Compute Bias Corrected Quantities.
#' @description Compute bias corrected residuals variance covariance matrix
#' and information matrix.
#' Also provides the leverage values and corrected sample size when adjust.n is set to TRUE.
#' @name estimate2
#' 
#' @keywords internal
.sscResiduals <- function(object, ssc, algorithm = "2"){
    algorithm <- match.arg(as.character(algorithm), choices = c("1","2"))

    ## ** initial values (i.e. non bias corrected)
    Omega0 <- ssc$Omega0 
    Omega.constrain <- ssc$Omega.constrain
    param0 <- ssc$param0
    residuals0 <- ssc$residuals0

    ## ** current values
    Omega <- object$sCorrect$moment$Omega 
    epsilon <- object$sCorrect$residuals
    leverage <- object$sCorrect$leverage
    dmu <- object$sCorrect$dmoment$dmu
    dOmega <- object$sCorrect$dmoment$dOmega
    vcov.param <- object$sCorrect$vcov.param
    
    endogenous <- object$sCorrect$endogenous
    n.endogenous <- length(endogenous)
    n.cluster <- object$sCorrect$cluster$n.cluster
    param.mean <- object$sCorrect$skeleton$Uparam.mean
    param.var <- object$sCorrect$skeleton$Uparam.var
    param.hybrid <- intersect(param.mean,param.var)
    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    n.pattern <- length(name.pattern)
    unique.pattern <- object$sCorrect$missing$unique.pattern

    if(length(param.mean)==0){
        stop("No mean parameter. No small sample correction needed. \n",
             "Consider setting  \'ssc\' to NA. \n")
    }
    ## ** Step (i-ii) compute individual and average bias
    dmu <- aperm(abind::abind(dmu[param.mean], along = 3), perm = c(3,2,1))
    vcov.muparam <- vcov.param[param.mean,param.mean,drop=FALSE]
    Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                  dimnames = list(endogenous, endogenous))
    n.Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                    dimnames = list(endogenous, endogenous))

    ## ls.Psi <- vector(mode = "list", length = n.cluster)
    for(iP in 1:n.pattern){ ## iP <- 1
        iY <- which(unique.pattern[iP,]>0)
        
        for(iC in missing.pattern[[iP]]){ ## iC <- 1
            ## individual bias
            if(length(param.mean)==1){
                iPsi <- vcov.muparam[1,1] * tcrossprod(dmu[,iY,iC])
            }else{
                iPsi <- t(dmu[,iY,iC])  %*% vcov.muparam %*% dmu[,iY,iC]
            }

            ## ls.Psi[[iC]] <- iPsi
            ## cumulated bias            
            Psi[iY,iY] <- Psi[iY,iY] + iPsi
            n.Psi[iY,iY] <- n.Psi[iY,iY] + 1
        }
    }
    ## take the average
    Psi[n.Psi>0] <- Psi[n.Psi>0]/n.Psi[n.Psi>0]
    
    ## ** Step (iii): compute corrected residuals and effective sample size
    ## done in estimate2 via moments2
    
    ## ** Step (iv): bias-corrected residual covariance matrix
    Omega.adj <- Omega0 + Psi

    ## ** Step (v): bias-corrected variance parameters
    A <- ssc$A

    ## *** right hand side of the equation
    index.upper.tri <- ssc$index.upper.tri[,"index"]
    eq.rhs <- stats::setNames((Omega.adj-Omega.constrain)[index.upper.tri],
                              ssc$name.rhs)

    ## *** left hand side of the equation
    index.Psi <- ssc$index.Psi

    if(NROW(index.Psi)>0){
        Z <- object$sCorrect$moment$iIB %*% object$sCorrect$moment$Lambda
        tZ <- t(Z)
        n.index.Psi <- NROW(index.Psi)
    
        ## A = t(Z) Psi Z + Sigma
        ## (t(Z) Psi Z)_{ij} = \sum_{k,l} Z_{k,i} Psi_{k,l} Z_{l,j}
        ## (t(Z) Psi Z)_{ij} regarding param_(k,l) = Z_{k,i} Z_{l,j}
        for(iPsi in 1:n.index.Psi){ ## iPsi <- 3
            iNamePsi <- rownames(index.Psi)[iPsi]
            iRowPsi <- index.Psi[iPsi,"row"]
            iColPsi <- index.Psi[iPsi,"col"]
            A[,iNamePsi] <- A[,iNamePsi] + (tZ[,index.Psi[iPsi,"row"]] %o% Z[index.Psi[iPsi,"col"],])[index.upper.tri]
            if(iRowPsi!=iColPsi){
                A[,iNamePsi] <- A[,iNamePsi] + (tZ[,index.Psi[iPsi,"col"]] %o% Z[index.Psi[iPsi,"row"],])[index.upper.tri]
            }
        }
        
    }

    ## *** solve equation
    ## microbenchmark::microbenchmark(svd = {asvd <- svd(A) ; asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs;},
    ## qr = qr.coef(qr(A), eq.rhs),
    ## Rcpp = OLS_cpp(A, eq.rhs),
    ## RcppTry = try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE),
    ## Rcpp2 = OLS2_cpp(A, eq.rhs),
    ## OLS1 = solve(crossprod(A), crossprod(A, eq.rhs)),
    ## OLS2 = solve(t(A) %*% A) %*% t(A) %*% eq.rhs,
    ## OLS_stats = stats::lsfit(x = A, y = eq.rhs),
    ## OLS_LINPACK = .Call(stats:::C_Cdqrls, x = A, y = eq.rhs, tolerance = 1e-7, FALSE)$coefficients, times = 500)
    if(lava.options()$method.estimate2=="svd"){
        asvd <- svd(A)
        iSolution <- try((asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs)[,1], silent = TRUE)
    }else if(lava.options()$method.estimate2=="ols"){
        iSolution <- try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE)
    }else{
        stop("unknown OLS methods \n")
    }
    
    if(inherits(iSolution, "try-error")){
        if(abs(det(t(A) %*% A)) <  1e-10){            
            stop("Singular matrix: cannot update the estimates \n")
        }else{
            stop(iSolution)
        }
    }
    names(iSolution) <- attr(A,"name")

    ## *** update parameters
    ## param0[ssc$name.var] - iSolution
    param0[names(iSolution)] <- iSolution

    ## ** Step (vi-vii): update derivatives and information matrix
    ## done in estimate2 via moments2
        
    ## ** Export
    attr(param0,"Omega") <- Omega.adj 
    attr(param0,"Psi") <- Psi
    return(param0)
}



##----------------------------------------------------------------------
### sCorrect-sscResiduals.R ends here

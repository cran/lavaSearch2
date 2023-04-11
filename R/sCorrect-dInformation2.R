### sCorrect-dInformation2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: Jan 17 2022 (18:44) 
##           By: Brice Ozenne
##     Update #: 349
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .dInformation2
#' @title Compute the First Derivative of the Expected Information Matrix
#' @description Compute the first derivative of the expected information matrix.
#' @name .dinformation2-internal
#' 
#' @details \code{calc_dinformation} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.dInformation2 <- function(dmu, dOmega, d2mu, d2Omega, OmegaM1,
                           missing.pattern, unique.pattern, name.pattern,
                           grid.3varD1, grid.2meanD1.1varD1, grid.2meanD2.1meanD1, grid.2varD2.1varD1,
                           name.param, leverage, n.cluster, weights){

    if(lava.options()$debug){cat(".dInformation2\n")}
    if(!is.null(weights)){
        stop(".dInformation2 does not support weights. \n")
    }
    symmetrize <- TRUE

    ## ** Prepare
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    n.grid.3varD1 <- NROW(grid.3varD1)
    if(symmetrize && n.grid.3varD1>0){
        grid.3varD1 <- grid.3varD1[grid.3varD1$duplicatedXYZ==FALSE,,drop=FALSE]
        n.grid.3varD1 <- NROW(grid.3varD1)
    }

    n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    if(symmetrize && n.grid.3varD1>0){
        grid.2meanD1.1varD1 <- grid.2meanD1.1varD1[grid.2meanD1.1varD1$duplicatedXY==FALSE,,drop=FALSE]
        n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    }

    n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)
    if(symmetrize && n.grid.2varD2.1varD1>0){
        grid.2varD2.1varD1 <- grid.2varD2.1varD1[grid.2varD2.1varD1$d2XZ+grid.2varD2.1varD1$d2YZ>0,,drop=FALSE]
        n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)

        grid.2varD2.1varD1[,"d2XZ"] <- grid.2varD2.1varD1[,"d2XZ"]*(1-grid.2varD2.1varD1[,"duplicatedXZ"])
        grid.2varD2.1varD1[,"d2YZ"] <- grid.2varD2.1varD1[,"d2YZ"]*(1-grid.2varD2.1varD1[,"duplicatedYZ"])
    }
    
    n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)
    if(symmetrize && n.grid.2meanD2.1meanD1>0){
        grid.2meanD2.1meanD1 <- grid.2meanD2.1meanD1[grid.2meanD2.1meanD1$d2XZ+grid.2meanD2.1meanD1$d2YZ>0,,drop=FALSE]
        n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)

        grid.2meanD2.1meanD1[,"d2XZ"] <- grid.2meanD2.1meanD1[,"d2XZ"]*(1-grid.2meanD2.1meanD1[,"duplicatedXZ"])
        grid.2meanD2.1meanD1[,"d2YZ"] <- grid.2meanD2.1meanD1[,"d2YZ"]*(1-grid.2meanD2.1meanD1[,"duplicatedYZ"])
    }

    dInfo <- array(0,
                   dim = c(n.param, n.param, n.param),
                   dimnames = list(name.param, name.param, name.param))
    
    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        if(!is.null(leverage)){
            iN.corrected <- length(iIndex) - colSums(leverage[iIndex,iY,drop=FALSE])
        }else{
            iN.corrected <- length(iIndex)
        }
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)

        ## *** 3 first derivative regarding the variance
        if(n.grid.3varD1>0){
            for(iGrid in 1:n.grid.3varD1){ # iGrid <- 1
                iName1 <- grid.3varD1[iGrid,"X"]
                iName2 <- grid.3varD1[iGrid,"Y"]
                iName3 <- grid.3varD1[iGrid,"Z"]

                ## term 1
                iDiag1 <- diag(iOmegaM1 %*% idOmega[[iName3]] %*% iOmegaM1 %*% idOmega[[iName2]] %*% iOmegaM1 %*% idOmega[[iName1]])
                iDiag2 <- diag(iOmegaM1 %*% idOmega[[iName2]] %*% iOmegaM1 %*% idOmega[[iName3]] %*% iOmegaM1 %*% idOmega[[iName1]])
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] - 1/2 * sum(iDiag1 * iN.corrected + iDiag2 * iN.corrected)

                ## symmetrize (XYZ = XZY = YXZ = YZX = ZXY = ZYX)
                if(symmetrize){
                    dInfo[iName1,iName3,iName2] <- dInfo[iName1,iName2,iName3]
                    dInfo[iName2,iName1,iName3] <- dInfo[iName1,iName2,iName3]
                    dInfo[iName2,iName3,iName1] <- dInfo[iName1,iName2,iName3]
                    dInfo[iName3,iName1,iName2] <- dInfo[iName1,iName2,iName3]
                    dInfo[iName3,iName2,iName1] <- dInfo[iName1,iName2,iName3]
                }
            }
        }

        ## *** 2 first derivative regarding the mean and one regarding the variance
        if(n.grid.2meanD1.1varD1>0){
            for(iGrid in 1:n.grid.2meanD1.1varD1){ # iGrid <- 1
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Y"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Z"]

                ## term 4
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] - sum(idmu[[iName1]] %*% iOmegaM1 %*% idOmega[[iName3]] %*% iOmegaM1 * idmu[[iName2]])
                
                ## symmetrize (XYZ = YXZ)
                if(symmetrize){
                    dInfo[iName2,iName1,iName3] <- dInfo[iName1,iName2,iName3]
                }

            }
        }
        
        ## *** 1 second derivative and 1 first derivative regarding the variance
        if(n.grid.2varD2.1varD1>0){
            for(iGrid in 1:n.grid.2varD2.1varD1){ # iGrid <- 1
                iName1 <- grid.2varD2.1varD1[iGrid,"X"]
                iName2 <- grid.2varD2.1varD1[iGrid,"Y"]
                iName3 <- grid.2varD2.1varD1[iGrid,"Z"]
               
                ## term 2
                if(grid.2varD2.1varD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var2"]
                    iDiag <- 1/2 * sum(diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName1]]) * iN.corrected)
                    dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + iDiag

                    ## symmetrize (XYZ = XZY)
                    if(symmetrize && (iName2 != iName3)){
                        dInfo[iName1,iName3,iName2] <- dInfo[iName1,iName3,iName2] + iDiag
                    }
                }

                ## term 3
                if(grid.2varD2.1varD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var2"]
                    iDiag <- 1/2 * sum(diag(iOmegaM1 %*% idOmega[[iName2]] %*% iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]]) * iN.corrected)
                    dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + iDiag

                    ## symmetrize (XYZ = ZYX)
                    if(symmetrize && (iName1 != iName3)){
                        dInfo[iName3,iName2,iName1] <- dInfo[iName3,iName2,iName1] + iDiag
                    }
                }

            }
        }

        ## *** 1 second derivative and 1 first derivative regarding the mean
        if(n.grid.2meanD2.1meanD1>0){
            for(iGrid in 1:n.grid.2meanD2.1meanD1){ # iGrid <- 1
                iName1 <- grid.2meanD2.1meanD1[iGrid,"X"]
                iName2 <- grid.2meanD2.1meanD1[iGrid,"Y"]
                iName3 <- grid.2meanD2.1meanD1[iGrid,"Z"]

                ## term 5
                if(grid.2meanD2.1meanD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var2"]
                    iDiag <- sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName2]])
                    dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + iDiag

                    ## symmetrize (XYZ=ZYX)
                    if(symmetrize && (iName1 != iName3)){
                        dInfo[iName3,iName2,iName1] <- dInfo[iName3,iName2,iName1] + iDiag
                    }
                }

                ## term 6
                if(grid.2meanD2.1meanD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var2"]
                    iDiag <- sum(idmu[[iName1]] %*% iOmegaM1 * id2mu[[d2.Var1]][[d2.Var2]])
                    dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + iDiag

                    ## symmetrize (XYZ = XZY)
                    if(symmetrize && (iName2 != iName3)){
                        dInfo[iName1,iName3,iName2] <- dInfo[iName1,iName3,iName2] + iDiag
                    }
                }
            }
        }
    }
    
    ## dInfo.bis <- .old_dInformation2(dmu = dmu, dOmega = dOmega, d2mu = d2mu, d2Omega = d2Omega, OmegaM1 = OmegaM1,
    ##                                 missing.pattern = missing.pattern, unique.pattern = unique.pattern, name.pattern = name.pattern, 
    ##                                 grid.dInformation = expand.grid(X = name.param, Y = name.param, Z = name.param, duplicated = FALSE, stringsAsFactors = FALSE),
    ##                                 name.param = name.param, name.param.dInformation = name.param,
    ##                                 leverage = leverage, n.cluster = n.cluster)
    ## dInfo.ter <- .old_dInformation2(dmu = dmu, dOmega = dOmega, d2mu = d2mu, d2Omega = d2Omega, OmegaM1 = OmegaM1,
    ##                                 missing.pattern = missing.pattern, unique.pattern = unique.pattern, name.pattern = name.pattern, 
    ##                                 grid.dInformation = expand.grid(X = name.param[2:3], Y = name.param[2:3], Z = name.param[2:3], duplicated = FALSE, stringsAsFactors = FALSE),
    ##                                 name.param = name.param[2:3], name.param.dInformation = name.param[2:3],
    ##                                 leverage = leverage, n.cluster = n.cluster)

    ## dInfo.bis[,,1]-t(dInfo.bis[,,1])
    ## dInfo.bis[,,2]-t(dInfo.bis[,,2])
    ## dInfo.bis[,,3]-t(dInfo.bis[,,3])
    ## print(range(dInfo.bis - dInfo))
    ## print(range(dInfo.ter  - dInfo))
    
    ## ** export
    return(dInfo)
}

## * .dVcov.param
.dVcov.param <- function(vcov.param, dInformation, n.param, name.param){
    dVcov.param <- array(0,
                         dim = c(n.param,n.param,n.param),
                         dimnames = list(name.param,name.param,name.param)
                         )
    
    for(iP in name.param){ ## iP <- "Y1"
        if(any(dInformation[,,iP]!=0)){
            dVcov.param[,,iP] <- - vcov.param %*% dInformation[,,iP] %*% vcov.param
        }
    }
    
    return(dVcov.param)
}

## * .dRvcov.param
.dRvcov.param <- function(score, hessian, vcov.param, dVcov.param, n.param, name.param){

    dRvcov.param <- array(0,
                          dim = c(n.param,n.param,n.param),
                          dimnames = list(name.param,name.param,name.param)
                          )
    
    score2_vcov.param <- crossprod(score) %*% vcov.param
    score_vcov.param <- score %*% vcov.param
        
    for(iP in name.param){ ## iP <- 1
        if(any(dVcov.param[,,iP]!=0)){
            term1 <- dVcov.param[,,iP] %*% score2_vcov.param
        }else{
            term1 <- matrix(0, nrow = n.param, ncol = n.param)
        }
        term2 <- vcov.param %*% hessian[iP,,] %*% score_vcov.param
        ## dRvcov.param[,,iP] <- term2 + t(term2) ## (what was lavaSearch2 doing before version 2.0.0)
        dRvcov.param[,,iP] <- term1 + t(term1) + term2 + t(term2)
    }

    return(dRvcov.param)
}



## * .old_dInformation2
## .old_dInformation2 <- function(dmu, dOmega, d2mu, d2Omega, OmegaM1,
##                                missing.pattern, unique.pattern, name.pattern, 
##                                grid.dInformation, name.param, name.param.dInformation,
##                                leverage, n.cluster){

##     if(lava.options()$debug){cat(".dInformation2\n")}

##     ## ** Prepare
##     n.grid <- NROW(grid.dInformation)
##     n.param <- length(name.param)
##     n.param.dInformation <- length(name.param.dInformation)
##     n.pattern <- length(name.pattern)

##     dInfo <-  array(0,
##                     dim = c(n.param, n.param, n.param.dInformation),
##                     dimnames = list(name.param, name.param, name.param.dInformation))

##     index.duplicated <- which(grid.dInformation$duplicated)
##     index.Nduplicated <- setdiff(1:n.grid, index.duplicated)
##     ## grid.dInformation[index.Nduplicated,,drop=FALSE]
    
##     ## ** loop over missing data pattern
##     for(iP in 1:n.pattern){ ## iP <- 1
##         iPattern <- name.pattern[iP]
##         iOmegaM1 <- OmegaM1[[iPattern]]
##         iIndex <- missing.pattern[[iPattern]]
##         iY <- which(unique.pattern[iP,]==1)

##         if(!is.null(leverage)){
##             iN.corrected <- length(iIndex) - colSums(leverage[iIndex,iY,drop=FALSE])
##         }else{
##             iN.corrected <- length(iIndex)
##         }
##         for(iGrid in index.Nduplicated){ # iGrid <- 1
##             iName1 <- grid.dInformation[iGrid,"X"]
##             iName2 <- grid.dInformation[iGrid,"Y"]
##             iNameD <- grid.dInformation[iGrid,"Z"]
##             ## cat("* ", iNameD," ",iName1,"",iName2,"\n")

##             ## *** identify relevant terms
##             test.Omega1 <- !is.null(dOmega[[iNameD]]) && !is.null(dOmega[[iName1]]) && !is.null(dOmega[[iName2]])
##             test.Omega2a <- !is.null(d2Omega[[iNameD]][[iName1]]) && !is.null(dOmega[[iName2]])
##             test.Omega2b <- !is.null(d2Omega[[iName1]][[iNameD]]) && !is.null(dOmega[[iName2]])
##             test.Omega3a <- !is.null(d2Omega[[iNameD]][[iName2]]) && !is.null(dOmega[[iName1]])
##             test.Omega3b <- !is.null(d2Omega[[iName2]][[iNameD]]) && !is.null(dOmega[[iName1]])

##             test.mu1a <- !is.null(d2mu[[iNameD]][[iName1]]) && !is.null(dmu[[iName2]])
##             test.mu1b <- !is.null(d2mu[[iName1]][[iNameD]]) && !is.null(dmu[[iName2]])
##             test.mu2a <- !is.null(d2mu[[iNameD]][[iName2nn]]) && !is.null(dmu[[iName1]])
##             test.mu2b <- !is.null(d2mu[[iName2]][[iNameD]]) && !is.null(dmu[[iName1]])
##             test.mu3 <- !is.null(dOmega[[iNameD]]) && !is.null(dmu[[iName1]]) && !is.null(dmu[[iName2]])

##             if((test.Omega1 + test.Omega2a + test.Omega2b + test.Omega3a + test.Omega3b + test.mu1a + test.mu1b + test.mu2a + test.mu2b + test.mu3) == 0){
##                 next
##             }
##             ## *** extract quantities for computations 
##             if(test.mu1a){
##                 d2mu.D1 <- d2mu[[iNameD]][[iName1]][iIndex,iY,drop=FALSE]
##             }else if(test.mu1b){
##                 d2mu.D1 <- d2mu[[iName1]][[iNameD]][iIndex,iY,drop=FALSE]
##             }
##             if(test.mu2a){
##                 d2mu.D2 <- d2mu[[iNameD]][[iName2]][iIndex,iY,drop=FALSE]
##             }else if(test.mu2b){
##                 d2mu.D2 <- d2mu[[iName2]][[iNameD]][iIndex,iY,drop=FALSE]
##             }
##             if(test.Omega2a){
##                 d2Omega.D1 <- d2Omega[[iNameD]][[iName1]][iY,iY,drop=FALSE]
##             }else if(test.Omega2b){
##                 d2Omega.D1 <- d2Omega[[iName1]][[iNameD]][iY,iY,drop=FALSE]
##             }
##             if(test.Omega3a){
##                 d2Omega.D2 <- d2Omega[[iNameD]][[iName2]][iY,iY,drop=FALSE]
##             }else{
##                 d2Omega.D2 <- d2Omega[[iName2]][[iNameD]][iY,iY,drop=FALSE]
##             }
                
##             ## *** pre-compute 
##             if(!is.null(dOmega[[iName1]])){
##                 iOmegaM1.dOmega.1 <- iOmegaM1 %*% dOmega[[iName1]][iY,iY,drop=FALSE]
##             }
##             if(!is.null(dOmega[[iName2]])){
##                 iOmegaM1.dOmega.2 <- iOmegaM1 %*% dOmega[[iName2]][iY,iY,drop=FALSE]
##             }                    
##             if(!is.null(dOmega[[iNameD]])){
##                 iOmegaM1.dOmega.D <- iOmegaM1 %*% dOmega[[iNameD]][iY,iY,drop=FALSE]
##             }
            
##             ## *** evaluate contributions to dInformation
##             if(test.Omega1){

##                 iOmegaM1.dOmega.1 %*% iOmegaM1.dOmega.D %*% iOmegaM1.dOmega.2
                
##                 iDiag1 <- diag(iOmegaM1.dOmega.D %*% iOmegaM1.dOmega.2 %*% iOmegaM1.dOmega.1)
##                 iDiag2 <- diag(iOmegaM1.dOmega.2 %*% iOmegaM1.dOmega.D %*% iOmegaM1.dOmega.1)
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum((iDiag1 + iDiag2) * iN.corrected)
##                 ## if(iName1!=iName2){
##                 ##     cat(iName1,"/",iName2,"/",iNameD,": ",iDiag1," ",iDiag2,"\n")
##                 ## }
##             }

##             if(test.Omega3a || test.Omega3b){
##                 iDiag <- diag(iOmegaM1 %*% d2Omega.D2 %*% iOmegaM1.dOmega.1)
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)                
##             }
            
##             if(test.Omega2a || test.Omega2b){
##                 iDiag <- diag(iOmegaM1.dOmega.2 %*% iOmegaM1 %*% d2Omega.D1)
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)
##             }

##             if(test.mu1a || test.mu1b){
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% iOmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
##             }

##             if(test.mu2a || test.mu2b){
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% iOmegaM1 * d2mu.D2)
##             }
                    
##             if(test.mu3){
##                 dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% iOmegaM1.dOmega.D %*% iOmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
##             }

##         }
##     }

##     ### ** export
##     return(dInfo)
## }

######################################################################
### sCorrect-dInformation2.R ends here

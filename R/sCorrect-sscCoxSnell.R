### biasCoxSnell.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (10:20) 
## Version: 
## Last-Updated: Jan 11 2022 (17:36) 
##           By: Brice Ozenne
##     Update #: 177
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .init_sscCoxSnell
.init_sscCoxSnell <- function(object,...){
    return(object,...)
}

## * .sscCoxSnell
.sscCoxSnell <- function(object, ssc){
    param <- ssc$param0
    name.param <- names(param)
    
    ## ** compute JJK
    object$sCorrect$skeleton$grid.2varD2.1varD1
    JJK <- .calcJJK(object)

    ## print(range(.old_calcJJK(object)-JJK))
    ## print(sum(.old_calcJJK(object))-sum(JJK))

    ## ** least squares
    Y <- (1/2) * sapply(name.param, function(iP){sum(JJK[,,iP] * object$sCorrect$vcov.param)})
    X <- object$sCorrect$information
    ## information(e.lvm2)

    ## dd <- as.data.frame(cbind(Y = Y,X))
    ## names(dd) <- c("Y",paste0("X",1:9))
    ## e.lm <- lm(Y ~ -1+X1+X2+X3+X4+X5+X6+X7+X8+X9, data = dd)
    
    e.lm <- stats::lm.fit(y = Y, x = X)
    newparam <- param - e.lm$coefficient
    ## print(e.lm$coefficient)
    
    ## ** export
    attr(newparam,"JJK") <- JJK
    attr(newparam,"lm") <- e.lm
    return(newparam)

}


## * .calcJJK
.calcJJK <- function(object){

    ## ** extract information
    dmu <- object$sCorrect$dmoment$dmu
    d2mu <- object$sCorrect$d2moment$d2mu
    dOmega <- object$sCorrect$dmoment$dOmega
    d2Omega <- object$sCorrect$d2moment$d2Omega

    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    unique.pattern <- object$sCorrect$missing$unique.pattern
    n.pattern <- length(name.pattern)
    OmegaM1 <- object$sCorrect$moment$OmegaM1.missing.pattern
    
    name.param <- names(object$sCorrect$param)
    n.param <- length(name.param)
    n.cluster <- object$sCorrect$cluster$n.cluster
    
    grid.2meanD1.1varD1 <- object$sCorrect$skeleton$grid.2meanD1.1varD1
    grid.2meanD2.1meanD1 <- object$sCorrect$skeleton$grid.2meanD2.1meanD1
    grid.2varD2.1varD1 <- object$sCorrect$skeleton$grid.2varD2.1varD1
    n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)
    n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)

    ## ** prepare output    
    JJK <-  array(0, dim = c(n.param,n.param,n.param),
                  dimnames = list(name.param,name.param,name.param))

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)

        ## *** 1 second derivative and 1 first derivative regarding the variance
        if(n.grid.2varD2.1varD1>0){
            for(iGrid in 1:n.grid.2varD2.1varD1){ # iGrid <- 1
                iName1 <- grid.2varD2.1varD1[iGrid,"X"]
                iName2 <- grid.2varD2.1varD1[iGrid,"Y"]
                iName3 <- grid.2varD2.1varD1[iGrid,"Z"]

                ## term 1
                if(grid.2varD2.1varD1[iGrid,"d2XY"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2XY.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2XY.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName3]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    ## cat("a: ",iName1," ",iName2," ",iName3,"\n")
                }

                ## term 2
                if(grid.2varD2.1varD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName2]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    ## cat("b: ",iName1," ",iName2," ",iName3,"\n")
                }

                ## term 3
                if(grid.2varD2.1varD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName1]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + 1/2 * sum(iDiag * n.cluster)
                    ## cat(iGrid,") c: ",iName1," ",iName2," ",iName3," = ", 1/2 * sum(iDiag * n.cluster),"\n")
                }
            }
        }

        ## *** 1 second derivative and 1 first derivative regarding the mean
        if(n.grid.2meanD2.1meanD1>0){
            for(iGrid in 1:n.grid.2meanD2.1meanD1){ # iGrid <- 1
                iName1 <- grid.2meanD2.1meanD1[iGrid,"X"]
                iName2 <- grid.2meanD2.1meanD1[iGrid,"Y"]
                iName3 <- grid.2meanD2.1meanD1[iGrid,"Z"]

                ## term 4
                if(grid.2meanD2.1meanD1[iGrid,"d2XY"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2XY.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2XY.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName3]])
                }

                ## term 5
                if(grid.2meanD2.1meanD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName2]])
                }

                ## term 6
                if(grid.2meanD2.1meanD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName1]])
                }
            }
        }
        
        ## *** 2 first derivative regarding the mean and one regarding the variance
        if(n.grid.2meanD1.1varD1>0){
            for(iGrid in 1:n.grid.2meanD1.1varD1){ # iGrid <- 1

                ## term 7
                iName1 <- grid.2meanD1.1varD1[iGrid,"Z"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Y"]
                value <- sum(idmu[[iName2]] %*% iOmegaM1 %*% idOmega[[iName1]] %*% iOmegaM1 * idmu[[iName3]])
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + value

                ## term 8 
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Z"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Y"]
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value
                
                ## term 9
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Y"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Z"]
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value
            }
        }
        

    }

    ## sum(abs(JJK)>0)
    return(JJK)
}


## * .old_calcJJK
.old_calcJJK <- function(object){

    ## ** extract information
    dmu <- object$sCorrect$dmoment$dmu
    d2mu <- object$sCorrect$d2moment$d2mu
    dOmega <- object$sCorrect$dmoment$dOmega
    d2Omega <- object$sCorrect$d2moment$d2Omega

    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    unique.pattern <- object$sCorrect$missing$unique.pattern
    n.pattern <- length(name.pattern)
    OmegaM1 <- object$sCorrect$moment$OmegaM1.missing.pattern
    
    name.param <- names(object$sCorrect$param)
    n.param <- length(name.param)
    n.cluster <- object$sCorrect$cluster$n.cluster
    
    grid.2meanD1.1varD1 <- object$sCorrect$skeleton$grid.2meanD1.1varD1
    grid.2meanD2.1meanD1 <- object$sCorrect$skeleton$grid.2meanD2.1meanD1
    grid.2varD2.1varD1 <- object$sCorrect$skeleton$grid.2varD2.1varD1
    n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)
    n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)

    ## ** prepare output    
    JJK <-  array(0, dim = c(n.param,n.param,n.param),
                  dimnames = list(name.param,name.param,name.param))

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)

        for(iParam1 in 1:n.param){ ## iParam1 <- 1
            for(iParam2 in 1:n.param){ ## iParam2 <- 1
                for(iParam3 in 1:n.param){ ## iParam3 <- 1

                    iName1 <- name.param[iParam1]
                    iName2 <- name.param[iParam2]
                    iName3 <- name.param[iParam3]
                    
                    ## *** 1 second derivative and 1 first derivative regarding the variance

                    ## term 1
                    if(!is.null(idOmega[[iName3]]) && !is.null(id2Omega[[iName1]][[iName2]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName1]][[iName2]] %*% iOmegaM1 %*% idOmega[[iName3]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    }else if(!is.null(idOmega[[iName3]]) && !is.null(id2Omega[[iName2]][[iName1]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName2]][[iName1]] %*% iOmegaM1 %*% idOmega[[iName3]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    }

                    ## term 2
                    if(!is.null(idOmega[[iName2]]) && !is.null(id2Omega[[iName1]][[iName3]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName1]][[iName3]] %*% iOmegaM1 %*% idOmega[[iName2]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    }else if(!is.null(idOmega[[iName2]]) && !is.null(id2Omega[[iName3]][[iName1]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName3]][[iName1]] %*% iOmegaM1 %*% idOmega[[iName2]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                    }

                    ## term 3
                    if(!is.null(idOmega[[iName1]]) && !is.null(id2Omega[[iName2]][[iName3]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName2]][[iName3]] %*% iOmegaM1 %*% idOmega[[iName1]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + 1/2 * sum(iDiag * n.cluster)
                        ## cat("c: ",iName1," ",iName2," ",iName3," = ", 1/2 * sum(iDiag * n.cluster),"\n")
                    }else if(!is.null(idOmega[[iName1]]) && !is.null(id2Omega[[iName3]][[iName2]])){
                        iDiag <- diag(iOmegaM1 %*% id2Omega[[iName3]][[iName2]] %*% iOmegaM1 %*% idOmega[[iName1]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + 1/2 * sum(iDiag * n.cluster)
                        ## cat("c: ",iName1," ",iName2," ",iName3," = ", 1/2 * sum(iDiag * n.cluster),"\n")
                    }
        
                    ## *** 1 second derivative and 1 first derivative regarding the mean

                    ## term 4
                    if(!is.null(idmu[[iName3]]) && !is.null(id2mu[[iName1]][[iName2]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[iName1]][[iName2]] %*% iOmegaM1 * idmu[[iName3]])
                    }else if(!is.null(idmu[[iName3]]) && !is.null(id2mu[[iName2]][[iName1]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[iName2]][[iName1]] %*% iOmegaM1 * idmu[[iName3]])
                    }

                    ## term 5
                    if(!is.null(idmu[[iName2]]) && !is.null(id2mu[[iName1]][[iName3]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[iName1]][[iName3]] %*% iOmegaM1 * idmu[[iName2]])
                    }else if(!is.null(idmu[[iName2]]) && !is.null(id2mu[[iName3]][[iName1]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[iName3]][[iName1]] %*% iOmegaM1 * idmu[[iName2]])
                    }

                    ## term 6
                    if(!is.null(idmu[[iName1]]) && !is.null(id2mu[[iName2]][[iName3]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + sum(id2mu[[iName2]][[iName3]] %*% iOmegaM1 * idmu[[iName1]])
                    }else if(!is.null(idmu[[iName1]]) && !is.null(id2mu[[iName3]][[iName2]])){
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + sum(id2mu[[iName3]][[iName2]] %*% iOmegaM1 * idmu[[iName1]])
                    }

                    ## *** 2 first derivative regarding the mean and one regarding the variance

                    ## term 7
                    if(!is.null(idmu[[iName2]]) && !is.null(idOmega[[iName1]]) && !is.null(idmu[[iName3]])){
                        value <- sum(idmu[[iName2]] %*% iOmegaM1 %*% idOmega[[iName1]] %*% iOmegaM1 * idmu[[iName3]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + value
                    }

                    ## term 8 
                    if(!is.null(idmu[[iName1]]) && !is.null(idOmega[[iName2]]) && !is.null(idmu[[iName3]])){
                        value <- sum(idmu[[iName1]] %*% iOmegaM1 %*% idOmega[[iName2]] %*% iOmegaM1 * idmu[[iName3]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value
                    }
                
                    ## term 9
                    if(!is.null(idmu[[iName2]]) && !is.null(idOmega[[iName3]]) && !is.null(idmu[[iName1]])){
                        value <- sum(idmu[[iName2]] %*% iOmegaM1 %*% idOmega[[iName3]] %*% iOmegaM1 * idmu[[iName1]])
                        JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value
                    }
                }
            }
        }        

    }

    return(JJK)
}
######################################################################
### biasCoxSnell.R ends here

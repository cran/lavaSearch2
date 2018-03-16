### skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (10:35) 
## Version: 
## Last-Updated: mar 12 2018 (17:59) 
##           By: Brice Ozenne
##     Update #: 798
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - skeleton
#' @title Pre-computation for the Score
#' @description Pre-compute quantities that are necessary to compute the score of a lvm model.
#' @name skeleton
#' 
#' @param object a \code{lvm} object.
#' @param df.param.all [data.frame] output of \code{\link{coefType}} containing the type of each coefficient.
#' @param param2originalLink [named character vector] matching between the name of the coefficient in lava and their label.
#' @param B,alpha.XGamma,Lambda,Psi [matrix] pre-computed matrix.
#' @param OD [list] the pre-computed quantities for the second derivatives. 
#' @param as.lava [logical] should the name of the links be used to name the coefficient?
#' Otherwise uses the labels (when defined) of each coefficient.
#' @param name.endogenous [character vector] name of the endogenous variables
#' @param name.latent [character vector] name of the latent variables
#' @param p [numeric vector, optional] vector of coefficients at which to evaluate the score.
#' @param data [data.frame, optional] data set.
#' @param ... [internal] only used by the generic method.
#' 
#' @details
#' When the use specify names for the coefficients (e.g. Y1[mu:sigma]) or uses constrains (Y1~beta*X1), \code{as.lava=FALSE} will use the names specified by the user (e.g. mu, sigma, beta) while \code{as.lava=TRUE} will use the name of the first link defining the coefficient.
#'
#' @examples
#' \dontrun{
#' skeleton <- lavaSearch2::skeleton
#' skeleton.lvm <- lavaSearch2::skeleton.lvm
#' skeleton.lvmfit <- lavaSearch2::skeleton.lvmfit
#' 
#' ## without constrain
#' m <- lvm(Y1~X1+X2+eta,Y2~X3+eta,Y3~eta)
#' latent(m) <- ~eta
#' 
#' e <- estimate(m,sim(m,1e2))
#' M.data <- as.matrix(model.frame(e))
#'
#' skeleton(e$model, as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), 
#'          update.value = FALSE)
#' skeleton(e, data = M.data, p = pars(e), as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), 
#'          update.value = TRUE)
#'
#' ## with constrains
#' m <- lvm(Y[mu:sigma] ~ beta*X1+X2)
#' e <- estimate(m, sim(m,1e2))
#' M.data <- as.matrix(model.frame(e))
#'
#' skeleton(e$model, as.lava = TRUE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, 
#'          update.value = FALSE)$skeleton
#' 
#' skeleton(e, data = M.data, p = pars(e), as.lava = FALSE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, 
#'          update.value = FALSE)$skeleton
#' 
#'}
#' @concept small sample inference
#' @concept derivative of the score equation
#' @keywords internal
`skeleton` <-
    function(object, ...) UseMethod("skeleton")


## * skeleton.lvm
#' @rdname skeleton
skeleton.lvm <- function(object, as.lava,
                         name.endogenous, name.latent,
                         ...){
    detail <- Y <- NULL ## [:for CRAN check] subset

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    
    ### ** prepare
    df.param.all  <- coefType(object, as.lava = FALSE)
    if(as.lava){
        param2originalLink <- subset(df.param.all, subset = !is.na(lava), select = c("originalLink", "param"))
        param2originalLink <- stats::setNames(param2originalLink$originalLink, param2originalLink$param)
    }else{
        param2originalLink <- subset(df.param.all, subset = !is.na(lava), select = "param", drop = TRUE)
        param2originalLink <- stats::setNames(param2originalLink, param2originalLink)
    }
    df.param.detail <- subset(df.param.all, subset = !is.na(detail))
    
    skeleton <- list()
    value <- list()

    ### ** Measurement model
    
    ## *** nu
    df.param.nu <-  subset(df.param.detail, subset = detail=="nu", select = c("value", "param", "Y", "name"))
    df.param.nu <- df.param.nu[order(df.param.nu$name),]
    value$nu <- stats::setNames(df.param.nu$value,df.param.nu$Y)
    skeleton$nu <- stats::setNames(param2originalLink[df.param.nu$param],df.param.nu$Y)

    ## *** X K
    df.param.K <- subset(df.param.detail, subset = detail == "K", select = c("value", "param", "X", "Y"))
    df.param.K <- df.param.K[order(df.param.K$Y),]

    if(NROW(df.param.K)>0){
        value$K <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){ # iEndogenous <- 1
            subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "value", drop = TRUE)
        }), name.endogenous)
            
        skeleton$K <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){
            param2originalLink[subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "param", drop = TRUE)]
        }), name.endogenous)
    
        skeleton$XK <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){
            subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "X", drop = TRUE)
        }), name.endogenous)
    }
    
    ## *** Lambda
    if(n.latent>0){
        ## define matrix
        value$Lambda <- matrix(0,nrow = n.latent, ncol = n.endogenous,
                               dimnames = list(name.latent,name.endogenous))
        skeleton$Lambda <- matrix(as.character(NA),nrow = n.latent, ncol = n.endogenous,
                                  dimnames = list(name.latent,name.endogenous))
        ## update according to the model
        df.param.Lambda <- subset(df.param.detail, subset = detail == "Lambda", select = c("X","Y","param","value","name"))
        df.param.Lambda$index <- match(df.param.Lambda$X, name.latent) + n.latent * (match(df.param.Lambda$Y, name.endogenous)-1)
        df.param.Lambda <- df.param.Lambda[order(df.param.Lambda$name),]

        ## store in the Lambda matrix the name of the coefficient and their pre-computed values
        dfNA.tempo <- subset(df.param.Lambda, subset = is.na(value))
        skeleton$Lambda[dfNA.tempo$index] <- stats::setNames(param2originalLink[dfNA.tempo$param],dfNA.tempo$Y)
        dfNNA.tempo <- subset(df.param.Lambda, subset = !is.na(value))
        value$Lambda[dfNNA.tempo$index] <- stats::setNames(dfNNA.tempo$value,dfNNA.tempo$Y)
        value$Lambda[!is.na(skeleton$Lambda)] <- NA
    }
    
    ## *** Sigma    
    ## define matrix
    value$Sigma <- matrix(0,nrow = n.endogenous, ncol = n.endogenous,
                          dimnames = list(name.endogenous,name.endogenous))
    skeleton$Sigma <- matrix(as.character(NA),nrow = n.endogenous, ncol = n.endogenous,
                             dimnames = list(name.endogenous,name.endogenous))

    ## update according to the model
    df.param.Sigma <- subset(df.param.detail,
                             subset = detail %in% c("Sigma_var","Sigma_cov"),
                             select = c("X","Y","param","value","name"))
    df.param.Sigma$index <- match(df.param.Sigma$X, name.endogenous) + n.endogenous*(match(df.param.Sigma$Y, name.endogenous)-1)

    dfNA.tempo <- subset(df.param.Sigma, subset = is.na(value))
    skeleton$Sigma[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
    dfNNA.tempo <- subset(df.param.Sigma, subset = !is.na(value))
    value$Sigma[dfNNA.tempo$index] <- dfNNA.tempo$value

    ## symmetrize
    skeleton$Sigma <- symmetrize(skeleton$Sigma, update.upper = TRUE)
    value$Sigma <- symmetrize(value$Sigma, update.upper = TRUE)
    value$Sigma[!is.na(skeleton$Sigma)] <- NA

### ** Structural model
    if(n.latent>0){
        ## *** alpha 
        df.param.alpha <-  subset(df.param.detail,
                                  subset = detail=="alpha",
                                  select = c("value","param","Y"))
        value$alpha <- stats::setNames(df.param.alpha$value,df.param.alpha$Y)
        skeleton$alpha <- param2originalLink[stats::setNames(df.param.alpha$param,df.param.alpha$Y)]

        ## *** X Gamma
        df.param.Gamma <- subset(df.param.detail,
                                 subset = detail=="Gamma",
                                 select = c("value","param","X","Y"))
        
        if(NROW(df.param.Gamma)>0){
            
            value$Gamma <- stats::setNames(lapply(1:n.latent, function(iLatent){
                subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "value", drop = TRUE)
            }), name.latent)
            
            skeleton$Gamma <- stats::setNames(lapply(1:n.latent, function(iLatent){ # iLatent <- 1
                param2originalLink[subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "param", drop = TRUE)]
            }), name.latent)
    
            skeleton$XGamma <- stats::setNames(lapply(1:n.latent, function(iLatent){
                subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "X", drop = TRUE)
            }), name.latent)
        }
        
        ## *** B
        ## define matrix
        value$B <- matrix(0,nrow = n.latent, ncol = n.latent,
                          dimnames = list(name.latent,name.latent))
        skeleton$B <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                             dimnames = list(name.latent,name.latent))

        if(any("B" %in% df.param.all$detail)){
            ## update according to the model
            df.param.B <- subset(df.param.detail,
                                 subset = detail == "B",
                                 select = c("X", "Y", "param", "value", "name"))
            df.param.B$index <- match(df.param.B$X, name.latent) + n.latent*(match(df.param.B$Y, name.latent)-1)
            dfNA.tempo <- subset(df.param.B, subset = is.na(value))            
            skeleton$B[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
            dfNNA.tempo <- subset(df.param.B, subset = is.na(value))
            value$B[dfNNA.tempo$index] <- dfNNA.tempo$value
            value$B[!is.na(skeleton$B)] <- NA            
        }
    
        ## *** Psi    
        ## define matrix
        value$Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                            dimnames = list(name.latent,name.latent))
        skeleton$Psi <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                               dimnames = list(name.latent,name.latent))
    
        ## update according to the model
        df.param.Psi <- subset(df.param.all,
                               subset = detail %in% c("Psi_var","Psi_cov"),
                               select = c("X", "Y", "param", "value", "Y", "name"))
                               
        df.param.Psi$index <- match(df.param.Psi$X, name.latent) + n.latent*(match(df.param.Psi$Y, name.latent)-1)

        dfNA.tempo <- subset(df.param.Psi, subset = is.na(value))      
        skeleton$Psi[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
        dfNNA.tempo <- subset(df.param.Psi, subset = is.na(value))
        value$Psi[dfNNA.tempo$index] <- dfNNA.tempo$value

        ## symmetrize
        skeleton$Psi <- symmetrize(skeleton$Psi, update.upper = TRUE)
        value$Psi <- symmetrize(value$Psi, update.upper = TRUE)
        value$Psi[!is.na(skeleton$Psi)] <- NA
    }

    ### ** export
    return(list(skeleton = skeleton,
                value = value,
                df.param = df.param.all,
                param2originalLink = param2originalLink)
           )
}


## * skeleton.lvmfit
#' @rdname skeleton
skeleton.lvmfit <- function(object, skeleton,
                            param, data,
                            name.endogenous, name.latent,
                            ...){
    ## note: object is not used except for finding the appropriate skeleton method.
    
    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    n.data <- NROW(data)

    ### ** Update skeleton with the current values

    ## *** nu
    index.update <- which(!is.na(skeleton$skeleton$nu))
    skeleton$value$nu[index.update] <- param[skeleton$skeleton$nu[index.update]]

    ## *** K 
    for(iY in 1:n.endogenous){ # iY <- 3
        if(length(skeleton$skeleton$K[[iY]])>0){
            index.update <- which(!is.na(skeleton$skeleton$K[[iY]]))
            skeleton$value$K[[iY]][index.update] <- param[skeleton$skeleton$K[[iY]][index.update]]
        }
    }

    ## *** Lambda
    if(n.latent>0){
        index.update <- which(!is.na(skeleton$skeleton$Lambda))
        skeleton$value$Lambda[index.update] <- param[skeleton$skeleton$Lambda[index.update]]
    }
    
    ## *** Sigma
    index.update <- which(!is.na(skeleton$skeleton$Sigma))
    skeleton$value$Sigma[index.update] <- param[skeleton$skeleton$Sigma[index.update]]

    ## *** linear predictor
    skeleton$value$nu.XK <- matrix(NA, nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                             dimnames = list(NULL,name.endogenous))
    for(iY in 1:n.endogenous){ # iY <- 1
        iY2 <- name.endogenous[iY]
        if(length(skeleton$value$K[[iY2]])>0){
            skeleton$value$nu.XK[,iY2] <- skeleton$value$nu[iY2] + data[,skeleton$skeleton$XK[[iY2]],drop=FALSE] %*% skeleton$value$K[[iY2]]
        }else{
            skeleton$value$nu.XK[,iY2] <- skeleton$value$nu[iY2]
        }
    }
        
    ### ** Structural model
    if(n.latent>0){
        ## *** alpha
        index.update <- which(!is.na(skeleton$skeleton$alpha))
        skeleton$value$alpha[index.update] <- param[skeleton$skeleton$alpha[index.update]]

        ## *** Gamma
        for(iLatent in 1:n.latent){
            if(length(skeleton$skeleton$Gamma[[iLatent]])>0){
                index.update <- which(!is.na(skeleton$skeleton$Gamma[[iLatent]]))
                skeleton$value$Gamma[[iLatent]][index.update] <- param[skeleton$skeleton$Gamma[[iLatent]][index.update]]
            }
        }
        
        ## *** B
        index.update <- which(!is.na(skeleton$skeleton$B))
        skeleton$value$B[index.update] <- param[skeleton$skeleton$B[index.update]]

        ## *** Psi
        index.update <- which(!is.na(skeleton$skeleton$Psi))
        skeleton$value$Psi[index.update] <- param[skeleton$skeleton$Psi[index.update]]

        ## *** linear predictor
        skeleton$value$alpha.XGamma <- matrix(NA,nrow = n.data, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL,name.latent))
        
        for(iLatent in 1:n.latent){
            iLatent2 <- name.latent[iLatent]
            if(length(skeleton$value$Gamma[[iLatent2]])>0){
                skeleton$value$alpha.XGamma[,iLatent2] <- skeleton$value$alpha[iLatent2] + data[,skeleton$skeleton$XGamma[[iLatent2]],drop=FALSE] %*% skeleton$value$Gamma[[iLatent2]]
            }else{
                skeleton$value$alpha.XGamma[,iLatent2] <- skeleton$value$alpha[iLatent2]
            }
        }
    }

  
### ** Export
    return(skeleton)
}


## * skeletonDtheta
#' @rdname skeleton
`skeletonDtheta` <-
    function(object, ...) UseMethod("skeletonDtheta")

## * skeletonDtheta.lvm
#' @rdname skeleton
skeletonDtheta.lvm <- function(object, data,
                               df.param.all, param2originalLink,
                               name.endogenous, name.latent, ...){

    factitious <- marginal <- param <- value <- X <- Y <- NULL ## [:for CRAN check] subset

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)

    df.param <- subset(df.param.all, subset = is.na(value) & marginal == FALSE & factitious == FALSE)
    Utype.by.detail <- tapply(df.param$detail, df.param$param, function(x){length(unique(x))})
    if(any(Utype.by.detail>1)){
        stop("cannot constrain two coefficients of different types to be equal \n")
    }
    name.param <- subset(df.param, subset = !duplicated(param), select = param, drop = TRUE)
    n.param <- length(name.param)

    name.originalLink <- as.character(param2originalLink)

### ** prepare
    n.data <- NROW(data)
    name.data <- colnames(data)
    
    mean.param <- c("nu","K","alpha","Gamma","Lambda","B")
    vcov.param <- c("Sigma_var","Sigma_cov","Psi_var","Psi_cov","Lambda","B")    
    dmu <- list()
    dOmega <- list()
    dLambda <- list()
    dB <- list()
    dPsi <- list()

    type <- stats::setNames(vector(mode = "character", n.param),name.originalLink)
    toUpdate <- stats::setNames(vector(mode = "logical", n.param),name.originalLink)
    
    ### ** Compute derivative or prepare for the derivative
    for(iName in name.param){ # iName <- name.param[1]

        iName2 <- as.character(param2originalLink[iName])
        type[iName2] <- unique(subset(df.param, subset = param == iName, select = "detail", drop = TRUE))
        iY <- subset(df.param, subset = param %in% iName, select = Y, drop = TRUE)
        iX <- subset(df.param, subset = param %in% iName, select = X, drop = TRUE)

        ## *** derivative regarding the mean        
        if(type[iName2] %in% mean.param){            
            if(type[iName2]=="nu"){
                dmu[[iName2]] <- matrix(as.numeric(name.endogenous %in% iY),
                                        nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(NULL, name.endogenous))
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="K"){
                dmu[[iName2]] <- matrix(0, nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(NULL, name.endogenous))
                for(Y.tempo in unique(iY)){                    
                    dmu[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="alpha"){
                dmu[[iName2]] <- matrix(as.numeric(name.latent %in% unique(iY)), nrow = n.data, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL, name.latent))                
                toUpdate[iName2] <- TRUE
            }else if(type[iName2]=="Gamma"){
                dmu[[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL, name.latent))
                for(Y.tempo in unique(iY)){ # Y.tempo <- "eta"
                    dmu[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- TRUE
            }
        }
        
        ## *** derivative regarding the residual variance covariance
        if(type[iName2] %in% vcov.param){
            
            if(type[iName2]=="Sigma_var"){
                dOmega[[iName2]] <- matrix(0,
                                           nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                           dimnames = list(name.endogenous, name.endogenous))
                dOmega[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="Sigma_cov"){
                dOmega[[iName2]] <- matrix(0,
                                           nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                           dimnames = list(name.endogenous, name.endogenous))
                dOmega[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                dOmega[[iName2]][match(iY, name.endogenous) + (match(iX, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }
            
        }        

        ## *** matrices
        if(type[iName2]=="Lambda"){            
            dLambda[[iName2]] <- matrix(0,
                                        nrow = n.latent, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(name.latent, name.endogenous))
            dLambda[[iName2]][match(iX, name.latent) + (match(iY, name.endogenous) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="B"){
            dB[[iName2]] <- matrix(0,
                                   nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                   dimnames = list(name.latent, name.latent))
            dB[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="Psi_var"){
            dPsi[[iName2]] <- matrix(0,
                                     nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                     dimnames = list(name.latent, name.latent))
            dPsi[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="Psi_cov"){
            dPsi[[iName2]] <- matrix(0,
                                     nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                     dimnames = list(name.latent, name.latent))
            dPsi[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            dPsi[[iName2]][match(iY, name.latent) + (match(iX, name.latent) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        } 
    }
### ** export
    return(list(
        dmu = dmu,
        dOmega = dOmega,
        dLambda = dLambda,
        dB = dB,
        dPsi = dPsi,
        type = type,
        toUpdate = toUpdate
    ))
}


## * skeletonDtheta.lvmfit
#' @rdname skeleton
skeletonDtheta.lvmfit <- function(object, dtheta,
                                  name.endogenous, name.latent,
                                  B, alpha.XGamma, Lambda, Psi,
                                  ...){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)

### ** Update partial derivatives
    if(any(dtheta$toUpdate)){
        type2update <- dtheta$type[dtheta$toUpdate]
        
        dtheta$iIB <- solve(diag(1,n.latent,n.latent)-B)
        dtheta$alpha.XGamma.iIB <- alpha.XGamma %*% dtheta$iIB
        dtheta$iIB.Lambda <-  dtheta$iIB %*% Lambda    
        dtheta$Psi.iIB <- Psi %*% dtheta$iIB
        dtheta$tLambda.tiIB.Psi.iIB <- t(dtheta$iIB.Lambda) %*% dtheta$Psi.iIB
        
        ## *** mean coefficients
        type.meanparam <- type2update[type2update %in% c("alpha","Lambda","Gamma","B")]
        n.meanparam <- length(type.meanparam)
        name.meanparam <- names(type.meanparam)
        
        if(n.meanparam>0){
            for(iP in 1:n.meanparam){ # iP <- 1
                iType <- type.meanparam[iP]
                iName <- name.meanparam[iP]
            
                if(iType == "alpha"){
                    dtheta$dmu[[iName]] <- dtheta$dmu[[iName]] %*% dtheta$iIB.Lambda
                }else if(iType == "Gamma"){
                    dtheta$dmu[[iName]] <- dtheta$dmu[[iName]] %*% dtheta$iIB.Lambda 
                }else if(iType == "Lambda"){
                    dtheta$dmu[[iName]] <- dtheta$alpha.XGamma.iIB %*% dtheta$dLambda[[iName]]
                }else if(iType == "B"){
                    dtheta$dmu[[iName]] <- dtheta$alpha.XGamma.iIB %*% dtheta$dB[[iName]] %*% dtheta$iIB.Lambda
                }

                colnames(dtheta$dmu[[iName]]) <- name.endogenous
            }
        }

        ## *** variance-covariance coefficients
        type.vcovparam <- type2update[type2update %in% c("Psi_var","Psi_cov","Lambda","B")]
        n.vcovparam <- length(type.vcovparam)
        name.vcovparam <- names(type.vcovparam)

        if(n.vcovparam>0){
            for(iP in 1:n.vcovparam){ # iP <- 1
                iType <- type.vcovparam[iP]
                iName <- name.vcovparam[iP]
        
                if(iType %in% "Psi_var"){
                    dtheta$dOmega[[iName]] <-  t(dtheta$iIB.Lambda) %*% dtheta$dPsi[[iName]] %*% dtheta$iIB.Lambda
                }else if(iType %in% "Psi_cov"){
                    dtheta$dOmega[[iName]] <-  t(dtheta$iIB.Lambda) %*% dtheta$dPsi[[iName]] %*% dtheta$iIB.Lambda
                }else if(iType == "Lambda"){
                    dtheta$dOmega[[iName]] <- dtheta$tLambda.tiIB.Psi.iIB %*% dtheta$dLambda[[iName]]
                    dtheta$dOmega[[iName]] <- dtheta$dOmega[[iName]] + t(dtheta$dOmega[[iName]])
                }else if(iType == "B"){
                    dtheta$dOmega[[iName]] <- dtheta$tLambda.tiIB.Psi.iIB %*% dtheta$dB[[iName]] %*% dtheta$iIB.Lambda
                    dtheta$dOmega[[iName]] <- dtheta$dOmega[[iName]] + t(dtheta$dOmega[[iName]])
                }

                colnames(dtheta$dOmega[[iName]]) <- name.endogenous
                rownames(dtheta$dOmega[[iName]]) <- name.endogenous
            }
        }
        
    }

### ** Export
    return(dtheta)

}

## * skeletonDtheta2
#' @rdname skeleton
`skeletonDtheta2` <-
    function(object, ...) UseMethod("skeletonDtheta2")

## * skeletonDtheta2.lvm
#' @rdname skeleton
skeletonDtheta2.lvm <- function(object, data, df.param.all,
                                param2originalLink, name.latent, ...){

    detail <- factitious <- marginal <- param <- value <- Y <- NULL ## [:for CRAN check] subset
    
    df.param <- subset(df.param.all, is.na(value) & marginal == FALSE & factitious == FALSE)
    dfred.param <- subset(df.param, subset = !duplicated(param))
    
    n.latent <- length(name.latent)
    n.data <- NROW(data)

### ** identify all combinations of coefficients with second derivative
    grid.mean <- list()

    grid.mean$alpha.B <- .combinationDF(dfred.param,
                                        detail1 = "alpha", name1 = "alpha",
                                        detail2 = "B", name2 = "B")

    grid.mean$alpha.Lambda <- .combinationDF(dfred.param,
                                             detail1 = "alpha", name1 = "alpha",
                                             detail2 = "Lambda", name2 = "Lambda")

    grid.mean$Gamma.B <- .combinationDF(dfred.param,
                                        detail1 = "Gamma", name1 = "Gamma",
                                        detail2 = "B", name2 = "B")

    grid.mean$Gamma.Lambda <- .combinationDF(dfred.param,
                                             detail1 = "Gamma", name1 = "Gamma",
                                             detail2 = "Lambda", name2 = "Lambda")
    
    grid.mean$Lambda.B <- .combinationDF(dfred.param,
                                        detail1 = "Lambda", name1 = "Lambda",
                                        detail2 = "B", name2 = "B")

    grid.mean$B.B <- .combinationDF(dfred.param,
                                    detail1 = "B", name1 = "B1",
                                    detail2 = "B", name2 = "B2")

    n.mean <- lapply(grid.mean, NROW)
    

    grid.vcov <- list()
    
    grid.vcov$Psi.Lambda <- .combinationDF(dfred.param,
                                           detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                           detail2 = "Lambda", name2 = "Lambda")

    grid.vcov$Psi.B <- .combinationDF(dfred.param,
                                      detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                      detail2 = "B", name2 = "B")

    grid.vcov$Lambda.B <- .combinationDF(dfred.param,
                                         detail1 = "Lambda", name1 = "Lambda",
                                         detail2 = "B", name2 = "B")

    grid.vcov$Lambda.Lambda <- .combinationDF(dfred.param,
                                              detail1 = "Lambda", name1 = "Lambda1",
                                              detail2 = "Lambda", name2 = "Lambda2")

    grid.vcov$B.B <- .combinationDF(dfred.param,
                                    detail1 = "B", name1 = "B1",
                                    detail2 = "B", name2 = "B2")
    
    n.vcov <- lapply(grid.vcov, NROW)
    
### ** convert back to lava names
    grid.mean <- lapply(grid.mean, function(x){ ## x <- grid.mean[[2]]
        if(length(x)>0){
            x[,1] <- param2originalLink[x[,1]]
            x[,2] <- param2originalLink[x[,2]]
        }
        return(x)
    })

    grid.vcov <- lapply(grid.vcov, function(x){ ## x <- grid.vcov[[2]]
        if(length(x)>0){
            x[,1] <- param2originalLink[x[,1]]
            x[,2] <- param2originalLink[x[,2]]
        }
        return(x)
    })

### ** prepare export
    if(any(unlist(n.mean)>0)){
        xx <- lapply(grid.mean, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, xx)
        name.tempo <- as.character(unique(collapseGrid[[1]]))
        d2mu <- lapply(name.tempo, function(x){
            iIndex <- which(collapseGrid[[1]]==x)
            v <- vector(mode = "list", length(iIndex))
            names(v) <- collapseGrid[[2]][iIndex]
            return(v)
        })
        names(d2mu) <- name.tempo
    }else{
        d2mu <- list()
    }
    
    if(any(unlist(n.vcov)>0)){
        xx <- lapply(grid.vcov, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, xx)
        name.tempo <- as.character(unique(collapseGrid[[1]]))
        d2Omega <- lapply(name.tempo, function(x){
            iIndex <- which(collapseGrid[[1]]==x)
            v <- vector(mode = "list", length(iIndex))
            names(v) <- collapseGrid[[2]][iIndex]
            return(v)
        })
        names(d2Omega) <- name.tempo
    }else{
        d2Omega <- list()
    }
    
    ## ** prepare alpha.B and alpha.Lambda
    if(any(df.param$detail == "alpha")){
        name.alpha <- subset(df.param, subset = !duplicated(param) & detail == "alpha", select = "param", drop = TRUE)
        ls.Malpha <- list()
        for(iName in name.alpha){ # iName <- name.alpha[1]

            iParam <- df.param[df.param$name == iName, "param"]
            iY <- subset(df.param, subset = param %in% iParam, select = Y, drop = TRUE)
            ls.Malpha[[iName]] <- matrix(as.numeric(name.latent %in% unique(iY)),
                                         nrow = n.data, ncol = n.latent, byrow = TRUE,
                                         dimnames = list(NULL, name.latent))
            
        }
    }
    if(n.mean$alpha.B>0){
        for(iP in 1:n.mean$alpha.B){ ## iP <- 1
            iName1 <- grid.mean$alpha.B[iP,"alpha"]
            iName2 <- grid.mean$alpha.B[iP,"B"]
            
            d2mu[[iName1]][[iName2]] <- ls.Malpha[[iName1]]
        }
    }
    if(n.mean$alpha.Lambda>0){
        for(iP in 1:n.mean$alpha.Lambda){ ## iP <- 1
            iName1 <- grid.mean$alpha.Lambda[iP,"alpha"]
            iName2 <- grid.mean$alpha.Lambda[iP,"Lambda"]
            
            d2mu[[iName1]][[iName2]] <- ls.Malpha[[iName1]]
        }
    }
    
    ## ** Store X for Gamma
    if(n.mean$Gamma.Lambda>0){
        for(iP in 1:n.mean$Gamma.Lambda){ ## iP <- 1
            iName1 <- grid.mean$Gamma.Lambda[iP,"Gamma"]
            iName2 <- grid.mean$Gamma.Lambda[iP,"Lambda"]

            iParam <- df.param[df.param$name == iName1, "param"]
            iX <- subset(df.param.all, subset = param %in% iParam, select = "X", drop = TRUE)
            iY <- subset(df.param.all, subset = param %in% iParam, select = "Y", drop = TRUE)

            d2mu[[iName1]][[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                               dimnames = list(NULL, name.latent))
            for(Y.tempo in unique(iY)){
                d2mu[[iName1]][[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
            }
        }
    }
    
    if(n.mean$Gamma.B>0){
        for(iP in 1:n.mean$Gamma.B){ ## iP <- 1
            iName1 <- grid.mean$Gamma.B[iP,"Gamma"]
            iName2 <- grid.mean$Gamma.B[iP,"B"]
            
            iParam <- df.param[df.param$name == iName1, "param"]
            iX <- subset(df.param.all, subset = param %in% iParam, select = "X", drop = TRUE)
            iY <- subset(df.param.all, subset = param %in% iParam, select = "Y", drop = TRUE)

            d2mu[[iName1]][[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                               dimnames = list(NULL, name.latent))
            for(Y.tempo in unique(iY)){
                d2mu[[iName1]][[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
            }
        }
    }

### ** export
    return(list(grid.mean = grid.mean,
                n.mean = n.mean,                
                grid.vcov = grid.vcov,
                n.vcov = n.vcov,
                d2mu = d2mu,
                d2Omega = d2Omega,
                toUpdate = any(unlist(c(n.mean,n.vcov))>0)
                ))
}

## * skeletonDtheta2.lvmfit
#' @rdname skeleton
skeletonDtheta2.lvmfit <- function(object, dtheta, d2theta,
                                   name.endogenous, name.latent,
                                   B, Lambda, Psi,
                                   ...){
    
    n.endogenous <- length(name.endogenous)

### ** second order partial derivatives
    if(any(d2theta$toUpdate)){
        
        ## *** mean coefficients        
        if(d2theta$n.mean$alpha.B>0){
            for(iP in 1:d2theta$n.mean$alpha.B){ # iP <- 1
                iName1 <- d2theta$grid.mean$alpha.B[iP,"alpha"]
                iName2 <- d2theta$grid.mean$alpha.B[iP,"B"]
                d2theta$d2mu[[iName1]][[iName2]] <- d2theta$d2mu[[iName1]][[iName2]] %*% dtheta$iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB.Lambda
            }
        }
        
        if(d2theta$n.mean$alpha.Lambda>0){
            for(iP in 1:d2theta$n.mean$alpha.Lambda){ # iP <- 1
                iName1 <- d2theta$grid.mean$alpha.Lambda[iP,"alpha"]
                iName2 <- d2theta$grid.mean$alpha.Lambda[iP,"Lambda"]

                d2theta$d2mu[[iName1]][[iName2]] <- d2theta$d2mu[[iName1]][[iName2]] %*% dtheta$iIB %*% dtheta$dLambda[[iName2]]
                
            }
        }

        if(d2theta$n.mean$Gamma.B>0){
            for(iP in 1:d2theta$n.mean$Gamma.B){ # iP <- 1
                iName1 <- d2theta$grid.mean$Gamma.B[iP,"Gamma"]
                iName2 <- d2theta$grid.mean$Gamma.B[iP,"B"]

                d2theta$d2mu[[iName1]][[iName2]] <- d2theta$d2mu[[iName1]][[iName2]] %*% dtheta$iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB.Lambda
            }
        }        

        if(d2theta$n.mean$Gamma.Lambda>0){
            for(iP in 1:d2theta$n.mean$Gamma.Lambda){ # iP <- 1
                iName1 <- d2theta$grid.mean$Gamma.Lambda[iP,"Gamma"]
                iName2 <- d2theta$grid.mean$Gamma.Lambda[iP,"Lambda"]                
                d2theta$d2mu[[iName1]][[iName2]] <- d2theta$d2mu[[iName1]][[iName2]] %*% dtheta$iIB %*% dtheta$dLambda[[iName2]]
            }
        }        

        if(d2theta$n.mean$Lambda.B>0){
            for(iP in 1:d2theta$n.mean$Lambda.B){ # iP <- 1
                iName1 <- d2theta$grid.mean$Lambda.B[iP,"Lambda"]
                iName2 <- d2theta$grid.mean$Lambda.B[iP,"B"]

                d2theta$d2mu[[iName1]][[iName2]] <- dtheta$alpha.XGamma.iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB %*% dtheta$dLambda[[iName1]]
            }
        }

        if(d2theta$n.mean$B.B>0){
            for(iP in 1:d2theta$n.mean$B.B){ # iP <- 1
                iName1 <- d2theta$grid.mean$B.B[iP,"B1"]
                iName2 <- d2theta$grid.mean$B.B[iP,"B2"]

                term1 <- dtheta$alpha.XGamma.iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB %*% dtheta$dB[[iName1]] %*% dtheta$iIB.Lambda
                term2 <- dtheta$alpha.XGamma.iIB %*% dtheta$dB[[iName1]] %*% dtheta$iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB.Lambda
                d2theta$d2mu[[iName1]][[iName2]] <- term1 + term2
            }
        }

        ## *** variance-covariance coefficients
        if(d2theta$n.vcov$Psi.Lambda>0){
            for(iP in 1:d2theta$n.vcov$Psi.Lambda){ # iP <- 1
                iName1 <- d2theta$grid.vcov$Psi.Lambda[iP,"Psi"]
                iName2 <- d2theta$grid.vcov$Psi.Lambda[iP,"Lambda"]

                term1 <- t(dtheta$dLambda[[iName2]]) %*% t(dtheta$iIB) %*% dtheta$dPsi[[iName1]] %*% dtheta$iIB.Lambda                
                d2theta$d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(d2theta$n.vcov$Psi.B>0){
            for(iP in 1:d2theta$n.vcov$Psi.B){ # iP <- 1
                iName1 <- d2theta$grid.vcov$Psi.B[iP,"Psi"]
                iName2 <- d2theta$grid.vcov$Psi.B[iP,"B"]

                term1 <- t(dtheta$iIB.Lambda) %*% t(dtheta$dB[[iName2]]) %*% t(dtheta$iIB) %*% dtheta$dPsi[[iName1]] %*% dtheta$iIB.Lambda
                d2theta$d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(d2theta$n.vcov$Lambda.B>0){
            for(iP in 1:d2theta$n.vcov$Lambda.B){ # iP <- 1
                iName1 <- d2theta$grid.vcov$Lambda.B[iP,"Lambda"]
                iName2 <- d2theta$grid.vcov$Lambda.B[iP,"B"]

                term1 <- t(dtheta$dLambda[[iName1]]) %*% t(dtheta$iIB) %*% t(dtheta$dB[[iName2]]) %*% t(dtheta$iIB) %*% Psi %*% dtheta$iIB.Lambda
                term2 <- t(dtheta$dLambda[[iName1]]) %*% t(dtheta$iIB) %*% Psi %*% dtheta$iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB.Lambda
                ## term2 <- dtheta$tLambda.tiIB.Psi.iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB %*% dtheta$dLambda[[iName1]]                
                d2theta$d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2)
            }
        }

        if(d2theta$n.vcov$Lambda.Lambda>0){
            for(iP in 1:d2theta$n.vcov$Lambda.Lambda){ # iP <- 1
                iName1 <- d2theta$grid.vcov$Lambda.Lambda[iP,"Lambda1"]
                iName2 <- d2theta$grid.vcov$Lambda.Lambda[iP,"Lambda2"]
                
                term1 <- t(dtheta$dLambda[[iName1]]) %*% t(dtheta$iIB) %*% dtheta$Psi.iIB %*% dtheta$dLambda[[iName2]]
                d2theta$d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(d2theta$n.vcov$B.B>0){
            for(iP in 1:d2theta$n.vcov$B.B){ # iP <- 1
                iName1 <- d2theta$grid.vcov$B.B[iP,"B1"]
                iName2 <- d2theta$grid.vcov$B.B[iP,"B2"]

                term1 <- t(dtheta$iIB.Lambda) %*% t(dtheta$dB[[iName2]]) %*% t(dtheta$iIB) %*% t(dtheta$dB[[iName1]]) %*% t(dtheta$iIB) %*% dtheta$Psi.iIB %*% Lambda
                term2 <- t(dtheta$iIB.Lambda) %*% t(dtheta$dB[[iName1]]) %*% t(dtheta$iIB) %*% t(dtheta$dB[[iName2]]) %*% t(dtheta$iIB) %*% dtheta$Psi.iIB %*% Lambda
                term3 <- t(dtheta$iIB.Lambda) %*% t(dtheta$dB[[iName1]]) %*% t(dtheta$iIB) %*% dtheta$Psi.iIB %*% dtheta$dB[[iName2]] %*% dtheta$iIB %*% Lambda
                d2theta$d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2) + term3 + t(term3)
            }
        }

    }

### ** Export
    return(d2theta)

}

## * .combination
#' @title Form all Unique Combinations Between two Vectors
#' @description Form all unique combinations between two vectors (removing symmetric combinations).
#' @name combination
#'
#' @param ... [vectors] elements to be combined.
#'
#' @return A matrix, each row being a different combination.
#' 
#' @examples
#' .combination <- lavaSearch2:::.combination
#' 
#' .combination(1,1)
#' .combination(1:2,1:2)
#' .combination(c(1:2,1:2),1:2)
#' 
#' .combination(alpha = 1:2, beta = 3:4)
#'
#' @keywords internal
.combination <- function(...){

    ## ** normalize arguments
    dots <- list(...)
    if(length(dots)!=2){
        stop("can only handle two vectors \n")
    }
    test.null <- unlist(lapply(dots,is.null))    
    if(any(test.null)){
        return(NULL)
    }
    dots <- lapply(dots,unique)

    ## ** form all combinations
    grid <- expand.grid(dots, stringsAsFactors = FALSE) 
    
    ## ** remove combinations (b,a) when (a,b) is already there
    name1 <- paste0(grid[,1],grid[,2])
    name2 <- paste0(grid[,2],grid[,1])

    if(NROW(grid)>1 && any(name1 %in% name2)){ 

        n.grid <- NROW(grid)
        test.duplicated <- c(FALSE,sapply(2:n.grid, function(iG){
            any(name2[iG] %in% name1[1:(iG-1)]) ## find duplicates
        }))

        grid <- grid[test.duplicated==FALSE,]
    }

    ## ** export
    return(grid)        
}


## * .combinationDF
.combinationDF <- function(data,
                           detail1, detail2,
                           name1, name2){

    detail <- NULL # [:for CRAN check] subset
    
    if(any(detail1 %in% data$detail) && any(detail2 %in% data$detail) ){
        ls.args <- list(subset(data, subset = detail %in% detail1, select = "param", drop = TRUE),
                        subset(data, subset = detail %in% detail2, select = "param", drop = TRUE))
        names(ls.args) <- c(name1,name2)
    
        return(do.call(.combination, args = ls.args))
        
    }else{
        
        return(numeric(0))
        
    }
}



##----------------------------------------------------------------------
### skeleton.R ends here


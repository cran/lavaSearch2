### coefType.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (14:38) 
## Version: 
## last-updated: jan 18 2018 (18:14) 
##           By: Brice Ozenne
##     Update #: 465
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - coefType
#' @title Extract the Specific Coefficient Bames or Positions in a LVM
#' @description Extract the specific coefficient names or positions in a latent varable model.
#' 
#' @name coefType
#' 
#' @param x a lvm model or a fitted lvm model 
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param data [optional] the dataset. Help to identify the categor
#' @param as.lava export the type of parameters mimicking lava:::coef.
#' @param keep.var should the variance parameters be returned?
#' @param ... arguments to be passed to \code{lava::coef}
#'
#' @details A lvm can be written as a measurement model:
#' \deqn{Y_i = \nu + \Lambda \eta_i + K X_i + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + B \eta_i + \Gamma X_i + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr
#'
#' \code{coefType} either returns the latin/greek letter corresponding to the coefficients
#' or it groups them:
#' \itemize{
#' \item intercept: \eqn{\nu} and \eqn{\alpha}.
#' \item regression: \eqn{\Lambda}, \eqn{K}, \eqn{B}, and \eqn{\Gamma}.
#' \item covariance: extra-diagonal terms of \eqn{\Sigma} and \eqn{\Psi}.
#' \item variance: diagonal of \eqn{\Sigma} and \eqn{\Psi}.
#' }
#'
#' @return \code{coefType} either returns a data.frame:
#' \itemize{
#' \item name: name of the link
#' \item Y: outcome variable
#' \item X: regression variable in the design matrix (could be a transformation of the original variables, e.g. dichotimization).
#' \item data: original variable
#' \item type: type of link
#' \item value: if TRUE, the value of the link is set and not estimated.
#' \item marginal: if TRUE, the value of the link does not impact the estimation.
#' \item detail: a more detailled decription of the type of link (see the details section)
#' \item lava: name of the parameter in lava
#' }
#' or a named vector containing the type of the links.
#' 
#' @examples 
#' #### regression ####
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE)
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#' 
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#'
#' coefReg(m)
#' coefReg(m, value = TRUE)
#' 
#' #### LVM ####
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1~y2
#'
#' m.Sim <- m
#' categorical(m.Sim, labels = c("a","b","c")) <- ~x2
#' e <- estimate(m, sim(m.Sim, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE) 
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#' 
#' coefExtra(m)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#' 
#' categorical(m, labels = as.character(1:3)) <- "X1"
#' coefExtra(m)
#' coefExtra(m, value = TRUE)
#'
#' categorical(m, labels = as.character(1:3)) <- "x1"
#'
#' coefType(m, as.lava = FALSE)
#' 
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#' coefIntercept(e)
#'
#' coefReg(e, value = TRUE)
#'
#' #### LVM with constrains ####
#' m <- lvm(c(Y1~0+1*eta1,Y2~0+1*eta1,Y3~0+1*eta1,
#'           Z1~0+1*eta2,Z2~0+1*eta2,Z3~0+1*eta2))
#' latent(m) <- ~eta1 + eta2
#' e <- estimate(m, sim(m,1e2))
#' coefType(m)
#' coef(e, level = 9)
#' coefType(e)
#' 
#' #### multigroup ####
#' m <- lvm(Y~X1+X2)
#' eG <- estimate(list(m,m), list(sim(m, 1e2),sim(m, 1e2)))
#' coefType(eG)
#' coef(eG)
#' 
#' coefIndexModel(eG)
#' 
#' @export
`coefType` <-
  function(x,...) UseMethod("coefType")


## * method - coefType
## ** method - coefType.lvm
#' @rdname coefType
#' @export
coefType.lvm <- function(x, data = NULL, as.lava = TRUE, ...){

    externalLink <- type <- NULL ## [:for CRAN check] subset
    
    ## *** extract all coef
    index.all <- which(!is.na(x$M), arr.ind = FALSE)
    ls.name <- list()
    ls.X <- list()
    ls.Y <- list()
    ls.type <- list()
    ls.value <- list()
    ls.param <- list()
    ls.marginal <- list()

    ## *** intercept
    n.intercept <- length(x$mean)
    if(n.intercept>0){
        ls.name$intercept <- names(x$mean)    
    
        ls.Y$intercept <- ls.name$intercept
        ls.X$intercept <- rep(NA, n.intercept)    
        ls.type$intercept <- rep("intercept", n.intercept)
        ls.value$intercept <- lapply(x$mean, function(iP){if(is.numeric(iP)){iP}else{NA}})
        ls.param$intercept <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                         iPar = unlist(x$mean),
                                         iFix = !is.na(ls.value$intercept),
                                         iName = ls.name$intercept)
                                     )
        ls.marginal$intercept <-  ls.name$intercept %in% exogenous(x)
    }
    
    ## *** regression
    arrIndex.regression <- which(x$M==1, arr.ind = TRUE)
    index.regression <- which(x$M==1, arr.ind = FALSE)
    n.regression <- length(index.regression)
    if(n.regression>0){
    
        ls.Y$regression <- colnames(x$M)[arrIndex.regression[,"col"]]
        ls.X$regression <- rownames(x$M)[arrIndex.regression[,"row"]]
        ls.name$regression <- paste0(ls.Y$regression,
                                     lava.options()$symbols[1],
                                     ls.X$regression)
    
        ls.type$regression <- rep("regression", n.regression)    
        ls.value$regression <- x$fix[index.regression]
        ls.param$regression <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                          iPar = x$par[index.regression],
                                          iFix = !is.na(ls.value$regression),
                                          iName = ls.name$regression)
                                      )
        ls.marginal$regression <- rep(FALSE,n.regression)
    }

    ## *** covariance
    M.cov <- x$cov
    M.cov[upper.tri(M.cov)] <- 0
    
    arrIndex.vcov <- which(M.cov==1, arr.ind = TRUE)
    index.vcov <- which(M.cov==1, arr.ind = FALSE)
    n.vcov <- length(index.vcov)
    if(n.vcov>0){
        
        Y.vcov <- colnames(x$cov)[arrIndex.vcov[,"col"]]
        X.vcov <- rownames(x$cov)[arrIndex.vcov[,"row"]]
    
        name.vcov <- paste0(Y.vcov,
                            lava.options()$symbols[2],
                            X.vcov)
    
        value.vcov <- x$covfix[index.vcov]
        param.vcov <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                 iPar = x$covpar[index.vcov],
                                 iFix = !is.na(value.vcov),
                                 iName = name.vcov)
                             )

        index.variance <- which(arrIndex.vcov[,1]==arrIndex.vcov[,2])
        ls.name$variance <- name.vcov[index.variance]
        n.variance <- length(ls.name$variance)
        ls.Y$variance <- Y.vcov[index.variance]
        ls.X$variance <- X.vcov[index.variance]
        ls.type$variance <- rep("variance", n.variance)
        ls.value$variance <- value.vcov[index.variance]
        ls.param$variance <- param.vcov[index.variance]
        ls.marginal$variance <- ls.name$variance %in% paste0(exogenous(x),lava.options()$symbols[2],exogenous(x))
    
        index.covariance <- which(arrIndex.vcov[,1]!=arrIndex.vcov[,2])
        ls.name$covariance <- name.vcov[index.covariance]
        n.covariance <- length(ls.name$covariance)
        ls.Y$covariance <- Y.vcov[index.covariance]
        ls.X$covariance <- X.vcov[index.covariance]
        ls.type$covariance <- rep("covariance", n.covariance)
        ls.value$covariance <- value.vcov[index.covariance]
        ls.param$covariance <- param.vcov[index.covariance]
        ls.marginal$covariance <- rep(FALSE, n.covariance)
    }
    
    ## *** external parameters
    n.external <- length(x$expar)
    if(n.external>0){
        ls.name$external <- names(x$expar)
        ls.type$external <- rep("external", n.external)

        ls.X$external <- rep(NA,n.external)
        for(iX in names(x$attributes$ordinalparname)){ ## iX <- "X1"
            ls.X$external[ls.name$external %in% x$attributes$ordinalparname[[iX]]] <- iX
        }
        ls.Y$external <- rep(NA,n.external)
        ls.value$external <- unlist(x$exfix)
        ls.param$external <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                        iPar = rep(NA,n.external),
                                        iFix = !is.na(ls.value$external),
                                        iName = ls.name$external)
                                    )
        ls.marginal$external <-  rep(FALSE, n.external)
    }

    ## *** merge
    df.param <- data.frame(name = unlist(ls.name),
                           Y = unlist(ls.Y),
                           X = unlist(ls.X),
                           data = unlist(ls.X),
                           type = unlist(ls.type),
                           value = unlist(ls.value),
                           param = unlist(ls.param),
                           marginal = unlist(ls.marginal),
                           stringsAsFactors = FALSE)
    df.param[df.param$X %in% latent(x),"data"] <- NA
    
    ## *** categorical variables
    if(!is.null(x$attributes$ordinalparname)){
        resCar <- defineCategoricalLink(x, link = df.param$name, data = data)
        
        ## normal parameters
        resCar.Nexternal <- subset(resCar,
                                   subset = is.na(externalLink),
                                   select = c("link","type","factice","level","originalLink","externalLink"))
        ## rename according to the main data frame
        match.tempo <- match(c("link","type"),
                             names(resCar.Nexternal))
        names(resCar.Nexternal)[match.tempo] <- c("name","distribution") 
        
        df.Nexternal <- merge(subset(df.param, subset = type != "external"),
                              resCar.Nexternal, by = "name")

        ## external parameters
        resCar.external <- subset(resCar,
                                  subset = !is.na(externalLink),
                                  select = c("link", "endogenous", "exogenous", "type", "factice", "level", "originalLink", "externalLink"))
        resCar.external$X <- paste0(resCar.external$exogenous,
                                    resCar.external$level)

        ## rename according to the main data frame
        match.tempo <- match(c("link","endogenous","exogenous","type"),
                             names(resCar.external))
        names(resCar.external)[match.tempo] <- c("name","Y","data","distribution")
        resCar.external$param <- resCar.external$name

        for(iCol in c("type","value","marginal")){ # iCol <- "type"
            name2col <- stats::setNames(df.param[[iCol]],df.param$name)
            resCar.external[,iCol] <- name2col[resCar.external$originalLink]
        }
        df.param <- rbind(resCar.external[,names(df.Nexternal),drop=FALSE],
                          df.Nexternal)
    }else{
        df.param$factice <- FALSE
        df.param$level <- as.character(NA)
        df.param$externalLink <-  as.character(NA)
        df.param$originalLink <- df.param$name
    }

    ## *** merge with lava
    coef.lava <- coef(x)
    name.coef <- names(coef.lava)

    index.keep <- which(df.param$type!="external" & df.param$factice == FALSE & df.param$marginal == FALSE)
    df.param$detail <- as.character(NA)
    df.param[index.keep, "detail"] <- detailName(x,
                                                 name.coef = df.param[index.keep, "name"],
                                                 type.coef = df.param[index.keep, "type"])
    df.param$lava <- name.coef[match(df.param$originalLink,coef.lava)]
    df.param <- df.param[order(df.param$type,df.param$detail,df.param$name),,drop=FALSE]
    rownames(df.param) <- NULL

    ## *** export
    if(as.lava){
        ## add extra parameter as links
        vec.extra <- unique(stats::na.omit(df.param$externalLink))
        if(length(vec.extra)>0){
            df.extra <- data.frame(name = vec.extra, type = "extra",
                                   lava = name.coef[match(vec.extra,coef.lava)],
                                   stringsAsFactors = FALSE)
            df.param <- rbind(df.param[,c("name", "type", "lava")],
                              df.extra)
        }
        
        ## 
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        out <- out[!duplicated(names(out))]
        return(out[coef.lava])    
    }else{
        return(df.param)
    }
}

## ** method - coefType.lvmfit
#' @rdname coefType
#' @export
coefType.lvmfit <- function(x, as.lava = TRUE, ...){ 

    ## *** find type of the coefficients in the original model
    df.param <- coefType(x$model0, as.lava = FALSE)
    
    ## *** export
    if(as.lava){
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        coef.lava <- names(stats::coef(x))
        return(out[coef.lava])    
    }else{
        return(df.param)
    }
}

## ** method - coefType.multigroup
#' @rdname coefType
#' @export
coefType.multigroup <- function(x, as.lava = TRUE, ...){

    n.model <- length(x$lvm)
    df.param <- NULL
    
    for(iModel in 1:n.model){ # iModel <- 2

        df.param <- rbind(df.param,
                          cbind(coefType(x$lvm[[iModel]], as.lava = FALSE), model = iModel)
                          )
      
    }

    df.param$name <- paste0(df.param$model,"@",df.param$name)
    
    ## *** export
    if(as.lava){        
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        return(out)
    }else{
        return(df.param)
    }
}

## * method - coefCov
#' @rdname coefType
#' @export
`coefCov` <-
  function(x,...) UseMethod("coefCov")

#' @rdname coefType
#' @export
coefCov.lvm <- function(x, value = FALSE, keep.var = FALSE, ...){

    res <- retainType(type = coefType(x, ...),
                      validType = c("covariance", if(keep.var){"variance"}else{NULL}),
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefCov.lvmfit <- coefCov.lvm

#' @rdname coefType
#' @export
coefCov.multigroup <- coefCov.lvm

## * method - coefExtra

#' @rdname coefType
#' @export
`coefExtra` <-
  function(x,...) UseMethod("coefExtra")

#' @rdname coefType
#' @export
coefExtra.lvm <- function(x, value = FALSE, ...){
    
    res <- retainType(type = coefType(x, ...),
                      validType = "extra",
                      value = value) 
    
    return(res)    
}

#' @rdname coefType
#' @export
coefExtra.lvmfit <- coefExtra.lvm

#' @rdname coefType
#' @export
coefExtra.multigroup <- coefExtra.lvm

## * method - coefIndexModel
#' @rdname coefType
#' @export
`coefIndexModel` <-
  function(x,...) UseMethod("coefIndexModel")

#' @rdname coefType
#' @export
coefIndexModel.lvm <- function(x, ...){
    name.coef <- stats::coef(x)
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

#' @rdname coefType
#' @export
coefIndexModel.lvmfit <- function(x, ...){
    name.coef <- names(stats::coef(x))
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

#' @rdname coefType
#' @export
coefIndexModel.multigroup <- function(x, ...){
    n.model <- length(x$lvm)

    ## new coef names
    allCoef <- x$name
    n.allCoef <- length(allCoef)
    index.AllCoef <- x$coef
  
    index <- stats::setNames(rep(NA, n.allCoef), allCoef)
  
    for(iModel in 1:n.model){ # iModel <- 1
        index[index.AllCoef[[iModel]]] <- iModel
    }
  
    #### export
    return(index)
}

#' @rdname coefType
#' @export
coefIndexModel.multigroupfit <- function(x, ...){
    out <- coefIndexModel(x$model)
    return(out[names(stats::coef(x))])
}
  
## * method - coefIntercept

#' @rdname coefType
#' @export
`coefIntercept` <-
  function(x,...) UseMethod("coefIntercept")

#' @rdname coefType
#' @export
coefIntercept.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "intercept",
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefIntercept.lvmfit <- coefIntercept.lvm

#' @rdname coefType
#' @export
coefIntercept.multigroup <- coefIntercept.lvm

## * method - coefRef
#' @rdname coefType
#' @export
`coefRef` <-
  function(x,...) UseMethod("coefRef")

#' @rdname coefType
#' @export
coefRef.lvmfit <- function(x, value = FALSE, ...){
    
    res <- retainType(type = attr(coefType(x, ...), "reference"),
                      validType = TRUE,
                      value = value)

    return(res)    
}

#' @rdname coefType
#' @export
`coefReg` <-
  function(x,...) UseMethod("coefReg")

#' @rdname coefType
#' @export
coefReg.lvm <- function(x, value = FALSE, ...){
    
     res <- retainType(type = coefType(x, ...),
                      validType = "regression",
                      value = value)

     return(res)
}

#' @rdname coefType
#' @export
coefReg.lvmfit <- coefReg.lvm

#' @rdname coefType
#' @export
coefReg.multigroup <- coefReg.lvm

## * method - coefVar

#' @rdname coefType
#' @export
`coefVar` <-
  function(x,...) UseMethod("coefVar")

#' @rdname coefType
#' @export
coefVar.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "variance",
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefVar.lvmfit <- coefVar.lvm

#' @rdname coefType
#' @export
coefVar.multigroup <- coefVar.lvm


#----------------------------------------------------------------------
### coefType.R ends here

## * retainType  (needed for coefCov/Latent/Ref)
retainType <- function(type, validType, value){
  index.var <- which(type %in% validType)
  
  if(length(index.var)>0){
      if(value){
          return(names(type)[index.var])
      }else{
          return(index.var)
      }
  }else{
      return(NULL)
  }
}

## * detailName (needed for coefType)
detailName <- function(x, name.coef, type.coef){
    ls.links <- initVarLinks(name.coef)
    index.loading <- setdiff(which(ls.links$var2 %in% latent(x)),
                             which(type.coef %in% c("covariance","variance")))
    if(length(index.loading)>0){
        type.coef[index.loading] <- "loading"
    }

    index.measurement <- which(ls.links$var1 %in% endogenous(x))
    if(length(index.measurement)>0){
        type.coef[index.measurement] <- as.character(factor(type.coef[index.measurement],
                                                       levels = c("intercept","regression","loading","covariance","variance"),
                                                       labels = c("nu","K","Lambda","Sigma_cov","Sigma_var")))
    }
    index.structural <- setdiff(1:length(type.coef),index.measurement)
    if(length(index.structural)>0){
        type.coef[index.structural] <- as.character(factor(type.coef[index.structural],
                                                      levels = c("intercept","regression","loading","covariance","variance"),
                                                      labels = c("alpha","Gamma","B","Psi_cov","Psi_var")))
    }
    return(type.coef)
}

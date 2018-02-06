### coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (17:29) 
## Version: 
## Last-Updated: feb  5 2018 (17:22) 
##           By: Brice Ozenne
##     Update #: 182
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .coef2
#' @title Export Mean and Variance Coefficients from a nlme Model
#' @description Export mean and variance coefficients from a nlme model.
#' @name coef2
#'
#' @param object a \code{gls} or \code{lme} object.
#'  
#' @details The variance coefficients that are exported are the residual variance of each outcome. 
#' This is \eqn{\sigma^2} for the first one and \eqn{k^2 \sigma^2} for the remaining ones.
#'
#' @return A numeric vector named with the names of the coefficient with three attributes:
#' \itemize{
#' \item mean.coef: the name of the mean coefficients.
#' \item var.coef: the name of the variance coefficients.
#' \item cor.coef:  the name of the correlation coefficients.
#' }
#'
#' @concept extractor
#' @keywords internal
`.coef2` <-
    function(object) UseMethod(".coef2")

## * .coef2.gls
#' @rdname coef2
.coef2.gls <- function(object){

     ## *** mean coefficients
    mean.coef <- stats::coef(object)

    ## *** variance coefficients
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(sigma2 = stats::sigma(object),
                      stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE)
                      )^2
    }else{
        var.coef <- c(sigma2 = stats::sigma(object)^2)
    }

    ## *** covariance coefficients
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }

    p <- c(mean.coef, cor.coef, var.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    return(p)
}




## * .coef2.lme
#' @rdname coef2
.coef2.lme <- function(object){

     ## *** mean coefficients
    mean.coef <- nlme::fixef(object)

    ## *** variance coefficients
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(sigma2 = stats::sigma(object),stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = stats::sigma(object)^2)
    }

    ## *** random effect coefficients
    random.coef <- as.double(nlme::getVarCov(object))    
    names(random.coef) <- paste0("ranCoef",1:length(random.coef))

     ## *** correlation coefficients
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }
    
    p <- c(mean.coef, cor.coef, var.coef, random.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    attr(p, "ran.coef") <- names(random.coef)
    return(p)
}

## * .getFormula2
`.getFormula2` <-
    function(object) UseMethod(".getFormula2")

## * .getFormula2.gls
.getFormula2.gls <- function(object){
    return(evalInParentEnv(object$call$model))
}

## * .getFormula2.lme
.getFormula2.lme <- function(object){
    return(evalInParentEnv(object$call$fixed))
}


## * .getGroups2
#' @title Reconstruct the Cluster Variable from a nlme Model
#' @description Reconstruct the cluster variable from a nlme model.
#' @name getGroup2
#'
#' @param object a \code{gls} or \code{lme} object.
#' @param cluster [vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param data [data.frame] the data set.
#' @param ... [internal] Only used by the generic method.
#'  
#' @details The variance coefficients that are exported are the residual variance of each outcome.
#' This is \eqn{\sigma^2} for the first one and \eqn{k^2 \sigma^2} for the remaining ones.
#'
#' @return A list containing:
#' \itemize{
#' \item cluster: the cluster index for each observation.
#' \item n.cluster: the number of clusters.
#' \item endogenous: to which endogenous variable each observation corresponds to.
#' \item n.endogenous: the number of endogenous variables.
#' \item index.obs: a vector to convert observations from the vector format to the matrix format.
#' }
#'
#' @concept extractor
#' @keywords internal
`.getGroups2` <-
    function(object, ...) UseMethod(".getGroups2")

## * .getGroups2.gls
#' @rdname getGroup2
.getGroups2.gls <- function(object, cluster, data, ...){

    ### ** get cluster
    cluster2 <- as.numeric(nlme::getGroups(object))
    if(length(cluster2)==0){ ## no correlation
        if(length(cluster) == 1 && is.character(cluster)){
            cluster2 <- as.numeric(as.factor(data[[cluster]]))
        }else{
            if(length(cluster)!=NROW(data)){
                stop("length of cluster and data do not match \n")
            }
            cluster2 <- as.numeric(as.factor(cluster))
        }
    }
    n.cluster <- length(unique(cluster2))

    ### ** get outcome
    if(!is.null(object$modelStruct$varStruct)){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else if(!is.null(object$modelStruct$corStruct)){        
        vec.rep0 <- unlist(attr(object$modelStruct$corStruct,"covariate"))
        if(is.factor(vec.rep0)){
            name.endogenous <- levels(vec.rep0)
        }else{
            name.endogenous <- unique(vec.rep0)
        }
        vec.rep <- vec.rep0[order(order(cluster2))]
    }else{
        vec.rep <- unlist(tapply(cluster2, cluster2, function(x){1:length(x)}))
        # vec.rep <- rep("1",NROW(data))
        name.endogenous <- unique(vec.rep)
    }
    vec.rep <- as.numeric(factor(vec.rep, levels = name.endogenous))
    n.endogenous <- length(name.endogenous)
    
    ### ** convert observations from the vector format to the matrix format
    index.obs <- cluster2+(vec.rep-1)*n.cluster

    ### ** export
    name.endogenous <- as.character(name.endogenous)

    return(list(cluster = cluster2,
                n.cluster = n.cluster,
                endogenous = vec.rep,
                name.endogenous = name.endogenous,
                n.endogenous = n.endogenous,
                index.obs = index.obs))
}

## * .getGroups2.lme
#' @name getGroup2
.getGroups2.lme <- function(object, ...){

    ## ** get cluster
    if(NCOL(object$groups)!=1){
        stop("cannot only handle one random effect \n")
    }
    cluster <- as.numeric(nlme::getGroups(object))
    n.cluster <- length(unique(cluster))

    ## ** get outcome
    if(!is.null(object$modelStruct$varStruct)){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else{
        if(!is.null(object$modelStruct$corStruct)){        
            vec.rep0 <- unlist(attr(object$modelStruct$corStruct,"covariate"))
        }else{
            vec.rep0 <- unlist(tapply(cluster,cluster,function(x){1:length(x)}))
        }
    
        if(is.factor(vec.rep0)){
            name.endogenous <- levels(vec.rep0)
        }else{
            name.endogenous <- unique(vec.rep0)
        }
        vec.rep <- vec.rep0[order(order(cluster))]
    }
    vec.rep <- as.numeric(factor(vec.rep, levels = name.endogenous))
    n.endogenous <- length(name.endogenous)
    
    ## ** convert observations from the vector format to the matrix format
    index.obs <- cluster+(vec.rep-1)*n.cluster

### ** export
    name.endogenous <- as.character(name.endogenous)
    
    return(list(cluster = cluster,
                n.cluster = n.cluster,
                endogenous = vec.rep,
                name.endogenous = name.endogenous,
                n.endogenous = n.endogenous,
                index.obs = index.obs))
}

## * .getVarCov2
#' @title Reconstruct the Marginal Variance Covariance Matrix from a nlme Model
#' @description Reconstruct the marginal variance covariance matrix from a nlme model.
#' @name getVarCov2
#'
#' @param object a \code{gls} or \code{lme} object
#' @param param [numeric vector] the mean and variance coefficients.
#' @param attr.param [character vector] the type of each coefficients (mean or variance).
#' @param endogenous [integer vector] the endogenous variable to which each observation corresponds.
#' @param name.endogenous [character vector] the name of the endogenous variables. 
#' @param n.endogenous [integer >0] the number of endogenous variables.
#' @param cluster [vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param n.cluster [integer >0] the number of groups.
#' @param ... [internal] Only used by the generic method.
#'  
#' @details The marginal variance covariance matrix for gls model is of the form:
#' 
#' \tabular{cccc}{
#' \eqn{\Sigma =} \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \sigma_2 \rho_{1,2}} \tab \eqn{\sigma^2 \sigma_3 \rho_{1,3}} \cr
#' \tab . \tab \eqn{\sigma^2 \sigma_2^2} \tab \eqn{\sigma^2 \sigma_3 \rho_{1,3}} \cr
#' \tab . \tab . \tab \eqn{\sigma^2 \sigma_3^2}
#' }
#'
#' The marginal variance covariance matrix for lme model is of the form:
#' @return A list containing:
#' \itemize{
#' \item Omega: the marginal variance covariance matrix for a full sample.
#' \item ls.indexOmega: a list containing for each sample the subset of endogenous variables available.
#' }
#' 
#' @concept extractor
#' @keywords internal
`.getVarCov2` <-
    function(object, ...) UseMethod(".getVarCov2")

## * .getVarCov2.gls
#' @rdname getVarCov2
.getVarCov2.gls <- function(object, param, attr.param,
                            endogenous, name.endogenous, n.endogenous,
                            cluster, n.cluster, ...){

    var.coef <- param[attr.param$var.coef]
    cor.coef <- param[attr.param$cor.coef]

    ## ** Diagonal terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if(length(name.other)>0){            
        sigma2.base <- stats::setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }else{
        sigma2.base <- stats::setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }
    template <- diag(as.double(sigma2.base),
                     nrow = n.endogenous, ncol = n.endogenous)
    
    ## ** Extra-diagonal terms
    if(length(cor.coef)>0){
        index.lower <- which(lower.tri(template))
        index.lower.arr <- which(lower.tri(template),arr.ind = TRUE)
        vec.sigma.tempo <- apply(index.lower.arr,1,function(x){prod(sqrt(sigma2.base[x]))})        
        template[index.lower] <- cor.coef*vec.sigma.tempo
        template <- symmetrize(template)
    }    
    dimnames(template) <- list(name.endogenous, name.endogenous)
    
    ## ** Index for each cluster    
    ls.indexOmega <- tapply(endogenous, cluster, function(iRep){iRep})

    ## ** export
    return(list(Omega = template,
                ls.indexOmega = ls.indexOmega))
}

## * .getVarCov2.lme
#' @rdname getVarCov2
.getVarCov2.lme <- function(object, param, attr.param,
                            endogenous, name.endogenous, n.endogenous,
                            cluster, n.cluster, ...){

    ## ** prepare with gls
    out <- .getVarCov2.gls(object, param = param, attr.param = attr.param,
                    endogenous = endogenous, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                    cluster = cluster, n.cluster = n.cluster)

    ## ** add contribution of the random effect
    ran.coef <- param[attr.param$ran.coef]
    out$Omega <- out$Omega + ran.coef

    ## ** export
    return(out)    
}




##----------------------------------------------------------------------
### coef2.R ends here

### effects2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  4 2019 (10:28) 
## Version: 
## Last-Updated: jan 24 2022 (12:03) 
##           By: Brice Ozenne
##     Update #: 384
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * effects2 (documentation)
#' @title Effects Through Pathways With Small Sample Correction 
#' @description Test whether a path in the latent variable model correspond to a null effect.
#' Similar to \code{lava::effects} but with small sample correction (if any).
#' So far it only work for a single path related two variable composed of one or two edges.
#' @name effects2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param linfct [character vector] The path for which the effect should be assessed (e.g. \code{"A~B"}),
#' i.e. the effect of the right variable (B) on the left variable (A). 
#' @param robust [logical] should robust standard errors be used instead of the model based standard errors? Should be \code{TRUE} if argument cluster is not \code{NULL}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param conf.level [numeric, 0-1] level of the confidence intervals.
#' @param from,to alternative to argument \code{linfct}. See \code{lava::effects}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/code{FALSE}/code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object.
#' 
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the confidence intervals.
#' 
#' @return A data.frame with a row per path.
#' 
#' @concept inference
#' @keywords smallSampleCorrection
#' @export
`effects2` <-
  function(object, linfct, robust, cluster, conf.level, ...) UseMethod("effects2")

## * effects2 (examples)
## TODO

## * effects2.lvmfit
#' @rdname effects2
#' @export
effects2.lvmfit <- function(object, linfct, robust = FALSE, cluster = NULL, conf.level = 0.95, to = NULL, from = NULL, df = lava.options()$df, ssc = lava.options()$ssc, ...){

    return(effects2(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...), linfct = linfct, to = to, from = from, robust = robust, cluster = cluster, conf.level = conf.level))

}

## * effects2.lvmfit2
#' @rdname effects2
#' @export
effects2.lvmfit2 <- function(object, linfct, robust = FALSE, cluster = NULL, conf.level = 0.95, to = NULL, from = NULL, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }
    object0 <- object
    class(object0) <- setdiff(class(object0),"lvmfit2")
        
    ## ** identify path
    if(!is.null(to) || !is.null(from)){
        n.hypo <- 1

        if(!missing(linfct)){
            stop("Cannot specify argument \'linfct\' at the same time as argument \'from\' or \'to\'. \n")
        }
        e.effects <- effects(object0, from = from, to = to)
        pathEffect <- stats::setNames(list(stats::setNames(list(e.effects$path), paste0(e.effects["to"],"~",e.effects["from"]))),paste0(e.effects["to"],"~",e.effects["from"]))
        type <- "total"
        null <- 0

    }else{
        n.hypo <- length(linfct)
        pathEffect <- vector(mode = "list", length = n.hypo)
        if(is.null(names(linfct))){
            names(pathEffect) <- linfct
        }else{
            if(any(duplicated(names(linfct)))){
                stop("Duplicated names for argument \'linfct\'. \n")
            }
            names(pathEffect) <- names(linfct)
        }
        type <- rep(as.character(NA), n.hypo)
        null <- rep(as.numeric(NA), n.hypo)

        for(iH in 1:n.hypo){

            if(grepl("|",linfct[iH], fixed=TRUE)){
                type[iH] <- base::trimws(strsplit(linfct[iH],split="|",fixed=TRUE)[[1]][2], which = "both")
                type[iH] <- match.arg(type[iH], c("indirect","direct","total"))
                linfct[iH] <- strsplit(linfct[iH],split="|",fixed=TRUE)[[1]][1]
            }else{
                type[iH] <- "total"
            }
            
            ## extract left and right side of the equation
            if(length(grep("=",linfct[iH]))>1){
                stop("Each element of argument \'linfct\' should contain at most one \'=\' sign.\n",
                     "Something like: coef1-2*coef2=0. \n")
            }
            iContrast <- createContrast(linfct[iH])
            if(iContrast$null==0){
                iContrast <- createContrast(linfct[iH], rowname.rhs = FALSE)
            }

            null[iH] <- unname(iContrast$null)
            iLHS.hypo_factor <- as.double(iContrast$contrast)
            iLHS.hypo_coef <- unname(colnames(iContrast$contrast))
            iN.param <- length(iLHS.hypo_coef)

            pathEffect[[iH]] <- stats::setNames(vector(mode = "list", length = iN.param), iLHS.hypo_coef)
            attr(pathEffect[[iH]], "factor") <- iLHS.hypo_factor
            
            for(iCoef in 1:iN.param){ ## iCoef <- 1
                ## extract all paths for each coefficient
                pathEffect[[iH]][[iCoef]] <- effects(object0, stats::as.formula(iLHS.hypo_coef[[iCoef]]))$path
                if(length(pathEffect[[iH]][[iCoef]])==0){
                    stop("Could not find path relative to coefficient ",iLHS.hypo_coef[[iCoef]]," (linfct=",linfct[iH],"). \n")
                }else if(type[iH]=="direct" && any(sapply(pathEffect[[iH]][[iCoef]],length)>2)){
                    pathEffect[[iH]][[iCoef]] <- pathEffect[[iH]][[iCoef]][sapply(pathEffect[[iH]][[iCoef]],length)==2]
                    if(length(pathEffect[[iH]][[iCoef]])==0){
                        stop("Could not find direct path relative to coefficient ",iLHS.hypo_coef[[iCoef]]," (linfct=",linfct[iH],"). \n")
                    }
                }else if(type[iH]=="indirect" && any(sapply(pathEffect[[iH]][[iCoef]],length)>2)){
                    pathEffect[[iH]][[iCoef]] <- pathEffect[[iH]][[iCoef]][sapply(pathEffect[[iH]][[iCoef]],length)>2]
                    if(length(pathEffect[[iH]][[iCoef]])==0){
                        stop("Could not find indirect path relative to coefficient ",iLHS.hypo_coef[[iCoef]]," (linfct=",linfct[iH],"). \n")
                    }
                }
            }
            
        }
    }

    ## ** extract information
    ## 0-order: param
    object.param <- coef(object, as.lava = FALSE)
    object.paramAll <- coef2(object, type = 9, labels = 1)[,"Estimate"]
    name.param <- names(object.param)
    n.param <- length(name.param)

    ## 1-order: score
    if(robust){
        object.score <- score(object, cluster = cluster, as.lava = FALSE)
    }

    ## 2-order: variance covariance
    object.vcov.param <- vcov(object, as.lava = FALSE)
    if(robust){
        object.rvcov.param <- vcov(object, robust = TRUE, cluster = cluster, as.lava = FALSE)
    }
    
    test.df <- (object$sCorrect$df == "satterthwaite")
    if(test.df){
        object.dVcov.param <- object$sCorrect$dVcov.param

        if(robust && (lava.options()$df.robust != 1)){

            if(!is.null(cluster)){ ## update derivative according to cluster
                object.dRvcov.param <- .dRvcov.param(score = object.score,
                                                     hessian = hessian2(object, cluster = cluster),
                                                     vcov.param = object.vcov.param,
                                                     dVcov.param = object.dVcov.param,
                                                     n.param = n.param,
                                                     name.param = name.param)
                                              
            }else{
                dRvcov.param <- object$sCorrect$dRvcov.param
            }
        }
    }

    coefEffect <- pathEffect

    ## ** identify coefficients corresponding to path
    for(iH in 1:n.hypo){ ## iH <- 1
        for(iCoef in 1:iN.param){ ## iCoef <- 1
            iN.path <- length(pathEffect[[iH]][[iCoef]])
            for(iPath in 1:iN.path){ ## iPath <- 1
                coefEffect[[iH]][[iCoef]][[iPath]] <- paste(pathEffect[[iH]][[iCoef]][[iPath]][-1], pathEffect[[iH]][[iCoef]][[iPath]][-length(pathEffect[[iH]][[iCoef]][[iPath]])], sep = lava.options()$symbols[1])
                if(any(coefEffect[[iH]][[iCoef]][[iPath]] %in% name.param) == FALSE){
                    stop("Incorrect path: ",paste(pathEffect[[iH]][[iCoef]][[iPath]], collapse="->"),"\n",
                         "Could not find coefficient: \"",paste(coefEffect[[iH]][[iCoef]][[iPath]][coefEffect[[iH]][[iCoef]][[iPath]] %in% name.param == FALSE], collapse = "\" \""),"\".\n")
                }
            }
        }
    }

    ## ** point estimate
    vec.beta <- stats::setNames(rep(NA, length = n.hypo), names(pathEffect))
    for(iH in 1:n.hypo){ ## iH <- 1
        iValue.param <- lapply(coefEffect[[iH]], function(iCoef){ ## for each coefficient (e.g. Y~E1 - Y~E2 = 0)
            iValue.path <- lapply(iCoef, function(iName){prod(object.paramAll[iName])}) ## get effect through each path corresponding to a coefficient (e.g. Y~E: Y~E and Y~X and X~E, i.e. \beta1 and \beta2*\beta3)
            return(do.call("sum", iValue.path)) ## return total effect (e.g. \beta1 + \beta2*\beta3)
        })
        if(is.null(attr(coefEffect[[iH]],"factor"))){
            vec.beta[iH] <- unlist(iValue.param)
        }else{
            vec.beta[iH] <- sum(attr(coefEffect[[iH]],"factor") * unlist(iValue.param))
        }
    }

    ## ** variance
    ## *** partial derivative
    dbeta.dtheta <- matrix(NA, nrow = n.hypo, ncol = n.param, dimnames = list(names(pathEffect), name.param))
    for(iH in 1:n.hypo){ ## iH <- 1
        iValue.param <- lapply(coefEffect[[iH]], function(iCoef){  ## iCoef <- coefEffect[[iH]][[1]] ## for each coefficient (e.g. Y~E1 - Y~E2 = 0)
            iDValue.path <- colSums(do.call(rbind,lapply(iCoef, function(iName){ ## iName <- iCoef[[1]] ## get derivative through each path corresponding to a coefficient (e.g. Y~E: Y~E and Y~X and X~E, i.e. \beta1 and \beta2*\beta3)
                iDeriv <- stats::setNames(rep(0, n.param), name.param)
                iDeriv[intersect(iName,name.param)] <- prod(object.paramAll[iName])/object.paramAll[intersect(iName,name.param)]
                return(iDeriv)
            })))
        })
        if(is.null(attr(coefEffect[[iH]],"factor"))){
            dbeta.dtheta[iH,] <- iValue.param[[1]]
        }else{
            dbeta.dtheta[iH,] <- attr(coefEffect[[iH]],"factor") %*% do.call(rbind,iValue.param)
        }   
    }

    if(robust){
        Mvcov.beta <- dbeta.dtheta %*% object.rvcov.param %*% t(dbeta.dtheta)
    }else{
        Mvcov.beta <- dbeta.dtheta %*% object.vcov.param %*% t(dbeta.dtheta)
    }

    ## ** compute df
    if(test.df){
        vec.df <- dfSigma(contrast = dbeta.dtheta,
                          score = object.score,
                          vcov = object.vcov.param,
                          rvcov = object.rvcov.param,
                          dVcov = object.dVcov.param,
                          dRvcov = object.dRvcov.param,
                          keep.param = dimnames(object.dVcov.param)[[3]],                            
                          type = if(robust){lava.options()$df.robust}else{1})
    }else{
        vec.df <- rep(0, n.hypo)
    }

    ## ** gather everything in glht object
    linfct2 <- diag(1, ncol = n.hypo, nrow = n.hypo)
    dimnames(linfct2) <- list(names(pathEffect),names(pathEffect))

    out <- list(model = object,
                linfct = linfct2,
                rhs = null,
                coef = vec.beta,
                vcov = Mvcov.beta,
                df = vec.df,
                alternative = "two.sided",
                type = NULL,
                robust = robust,
                ssc = object$sCorrect$ssc$type,
                grad = dbeta.dtheta,
                path = pathEffect,
                global = NULL)
    class(out) <- c("glht2","glht")

    ## ** export
    return(out)

}

## * effects.lvmfit2
#' @rdname effects2
#' @export
effects.lvmfit2 <- effects2.lvmfit2

######################################################################
### effects2.R ends here

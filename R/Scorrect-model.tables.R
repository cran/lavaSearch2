### Scorrect-model.tables.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan 12 2022 (09:41) 
## Version: 
## Last-Updated: Jan 12 2022 (11:45) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Estimate, Confidence Intervals, and P-value With Small Sample Correction
#' @description Extract estimate, standard error, confidence intervals and p-values associated to each coefficient of a latent variable model.
#' Similar to \code{lava::confint} but with small sample correction.
#' @name confint2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param robust [logical] should robust standard errors be used instead of the model based standard errors? Should be \code{TRUE} if argument cluster is not \code{NULL}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] when \code{TRUE} uses the same names as when using \code{stats::coef}.
#' @param transform [function] transformation to be applied.
#' @param conf.level [numeric, 0-1] level of the confidence intervals.
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
#' @return A data.frame with a row per coefficient.
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`model.tables2` <- function(object, robust, cluster, transform,
                            as.lava, conf.level, ...) UseMethod("model.tables2")

## * model.tables.lvmfit
#' @export
model.tables.lvmfit <- function(x, ...){

    x.confint <- confint(x, ...)
    out <- cbind(summary(x, ...)$coef, upper = NA, lower = NA)
    colnames(out)[1:4] <- c("estimate","se","statistic","p.value")
    out <- out[,c("estimate","se","lower","upper","statistic","p.value")]
    out[rownames(x.confint),"lower"] <- x.confint[,1]
    out[rownames(x.confint),"upper"] <- x.confint[,2]

    return(out)
}

## * model.tables2.lvmfit
#' @export
model.tables2.lvmfit <- function(object, robust = FALSE, cluster = NULL,
                                 transform = NULL, as.lava = TRUE, conf.level = 0.95, 
                                 ssc = lava.options()$ssc, df = lava.options()$df, ...){

    return(model.tables(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...),
                        robust = robust, cluster = cluster, as.lava = as.lava, conf.level = conf.level,
                        transform = transform))

}

## * model.tables.lvmfit2
#' @export
model.tables2.lvmfit2 <- function(object, robust = FALSE, cluster = NULL,
                            transform = NULL, as.lava = TRUE, conf.level = 0.95,  ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    ## ** new model parameters
    param <- coef(object, as.lava = as.lava)
    name.param <- names(param)
    n.param <- length(name.param)

    ## ** new Wald test
    type <- object$sCorrect$skeleton$type
    type <- type[!is.na(type$lava),]
    null <- stats::setNames(rep(0, n.param),name.param)
    if(any(type$detail %in% c("Sigma_var","Psi_var"))){
        param.var <- type[type$detail %in% c("Sigma_var","Psi_var"),"param"]
        if(as.lava){
            null[names(object$sCorrect$skeleton$originalLink2param)[match(param.var,object$sCorrect$skeleton$originalLink2param)]] <- NA
        }else{
            null[object$sCorrect$skeleton$originalLink2param[match(param.var,object$sCorrect$skeleton$originalLink2param)]] <- NA
        }
    }
    table.all <- compare2(object,
                          linfct = name.param,
                          rhs = null,
                          robust = robust,
                          cluster = NULL,
                          F.test = FALSE,
                          as.lava = FALSE,
                          sep = c("",""))

    tableS.all <- summary(table.all, test = multcomp::adjusted("none"), transform = transform, conf.level = conf.level, rowname.rhs = FALSE)$table2

    if(as.lava){
        tableS.all <- tableS.all[names(object$sCorrect$skeleton$originalLink2param),,drop=FALSE]
    }else{
        tableS.all <- tableS.all[as.character(object$sCorrect$skeleton$originalLink2param),,drop=FALSE]
    }
    return(tableS.all)

}

## *  model.tables.lvmfit2
#' @export
model.tables.lvmfit2 <- function(x, ...){ ## necessary as model.tables must have x has first argument
    return(model.tables2(x, ...))
}

##----------------------------------------------------------------------
### Scorrect-model.tables.R ends here

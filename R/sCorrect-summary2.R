
### sCorrect-summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: jan 23 2024 (10:26) 
##           By: Brice Ozenne
##     Update #: 554
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - summary2
#' @title Latent Variable Model Summary After Small Sample Correction
#' @description Summarize a fitted latent variable model.
#' Similar to \code{stats::summary} with small sample correction.
#' @name summary2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param digit [integer > 0] the number of decimal places to use when displaying the summary.
#' @param robust [logical] should robust standard errors be used instead of the model based standard errors? Should be \code{TRUE} if argument cluster is not \code{NULL}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/\code{FALSE}/\code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... [logical] arguments passed to lower level methods.
#' 
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#'
#' @details \code{summary2} is the same as \code{summary}
#' except that it first computes the small sample correction (but does not store it).
#' So if \code{summary2} is to be called several times,
#' it is more efficient to pre-compute the quantities for the small sample correction
#' using \code{sCorrect} and then call \code{summary2}.
#'
#' \code{summary2} returns an object with an element \code{table2} containing the estimates, standard errors, degrees of freedom,
#' upper and lower limits of the confidence intervals, test statistics, and p-values.
#' 
#' @examples
#' #### simulate data ####
#' m <- lvm(Y~X1+X2)
#' set.seed(10)
#' d <- lava::sim(m, 2e1)
#'
#' #### latent variable models ####
#' e.lvm <- estimate(m, data = d)
#' summary(e.lvm)$coef
#' 
#' summary2(e.lvm)
#' summary2(e.lvm, ssc = "none")
#' 
#' @concept small sample inference
#' @export
`summary2` <-
  function(object, robust, cluster, digit, ...) UseMethod("summary2")

## * summary2.lvmfit
#' @rdname summary2
#' @export
summary2.lvmfit <- function(object, robust = FALSE, cluster = NULL, digit = max(5, getOption("digit")), ssc = lava.options()$ssc, df = lava.options()$df, ...){

    return(summary(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...), robust = robust, cluster = NULL, digit = digit))

}

## * summary2.lvmfit2
#' @rdname summary2
#' @export
summary2.lvmfit2 <- function(object, robust = FALSE, cluster = NULL, digit = max(5, getOption("digit")), ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    ## ** table with se, df, confint, p-value for the corrected parameters
    tableS.all <- model.tables(object, robust = robust, cluster = cluster, as.lava = TRUE)
    name.param <- rownames(tableS.all)
    n.param <- length(name.param)

    ## ** get and  normalize lava summary
    object0 <- object
    class(object0) <- setdiff(class(object0),c("lvmfit2"))
    object.summary <- summary(object0, digits = digit)
    
    previous.summary <- object.summary$coef
    object.summary$coef <- tableS.all[name.param,c("estimate","se","statistic","df","p.value"),drop=FALSE]


    ## find digit
    vec.char <- setdiff(object.summary$coefmat[,"Estimate"],"")
    digit <- max(c(nchar(gsub(".","",vec.char,fixed = TRUE)))-1,1)

    ## ** update summary
    ## *** vcov
    object.summary$vcov <- attr(object$dVcov, "vcov.param")[name.param,name.param]    

    ## *** coef
    lava.rownames <- rownames(previous.summary)
    
    ## *** coefmat
    name.label0 <- trimws(rownames(CoefMat(object0, labels = 0, level = 9)), which = "both")
    index.titleVariance <- which(name.label0=="Residual Variances:")
    if(length(index.titleVariance)>0){
        ## rename variance parameters from Y to Y~~Y
        index.vcov <- (index.titleVariance+1):length(name.label0)
        index.var <- setdiff(index.vcov,grep("~~",name.label0,fixed=TRUE)) ## exclude covariance parameters that are already correctly named
        name.label0[index.var] <- paste0(name.label0[index.var],lava.options()$symbols[2],name.label0[index.var])
    }
    table.coefmat <- object.summary$coefmat
    if(object$sCorrect$df=="satterthwaite"){
        colnames(table.coefmat)[3:5] <- c("t-value","P-value","df")
    }else{
        colnames(table.coefmat)[3:5] <- c("Z-value","P-value","df")
    }

    ## mimic lava:::CoefMat (called by lava:::summary.lvmfit)
    table.coef <- object.summary$coef
    e2add <- format(round(table.coef[,"estimate"], max(1, digit - 1)), digits = digit - 1)
    e2add <- gsub(" NA","",e2add)
    sd2add <- format(round(table.coef[,"se"], max(1, digit - 1)), digits = digit - 1)
    sd2add <- gsub(" NA","",sd2add)
    df2add <- as.character(round(table.coef[,"df"],2))    
    df2add[is.na(df2add)] <- ""
    t2add <- format(round(table.coef[,"statistic"], max(1, digit - 1)), digits = digit - 1)
    t2add <- gsub(" NA","",t2add)

    p2add <- formatC(table.coef[,"p.value"], digits = digit - 1, format = "g",  preserve.width = "common", flag = "")
    p2add <- gsub(" NA","",p2add)
    p2add[table.coef[,4] < 1e-12] <- "  <1e-12"

    M2add <- cbind(e2add,sd2add,t2add,p2add,df2add)
    table.coefmat[,"df"] <- ""
    table.coefmat[match(rownames(table.coef), name.label0),] <- M2add

    table.coefmat[object.summary$coefmat[,4]=="",4] <- ""
    object.summary$coefmat <- table.coefmat

    ## ** Export
    if(robust){
        colnames(object.summary$coefmat)[2] <- "robust SE"
        colnames(object.summary$coef)[2] <- "robust SE"
    }

    ## ** gather all results in one table
    object.summary$table2 <- data.frame(matrix(NA, nrow = n.param, ncol = 7,
                                               dimnames = list(name.param,
                                                               c("estimate","se","df","lower","upper","statistic","p.value"))
                                               ), stringsAsFactors = FALSE)

    object.summary$table2$estimate <- tableS.all[name.param,"estimate"]
    object.summary$table2$se <- tableS.all[name.param,"se"]
    object.summary$table2$df <- tableS.all[name.param,"df"]
    object.summary$table2$lower <- object.summary$table2$estimate + object.summary$table2$se * stats::qt(p=0.025, df = object.summary$table2$df)
    object.summary$table2$upper <- object.summary$table2$estimate + object.summary$table2$se * stats::qt(p=0.975, df = object.summary$table2$df)
    object.summary$table2$statistic <- tableS.all[name.param,"statistic"]
    object.summary$table2$p.value <- tableS.all[name.param,"p.value"]

    ## ** export
    return(object.summary)
    
}

## * summary.lvmfit2
#' @rdname summary2
#' @export
summary.lvmfit2 <- summary2.lvmfit2


##----------------------------------------------------------------------
### Scorrect-summary2.R ends here


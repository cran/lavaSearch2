### summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: feb  6 2018 (16:45) 
##           By: Brice Ozenne
##     Update #: 180
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - summary
#' @title  Summary with Small Sample Correction
#' @description Summary with small sample correction.
#' @name summary
#'
#' @param object a \code{gls}, \code{lme} or \code{lvm} object.
#' @param digit [integer > 0] the number of digit to keep when displaying the summary.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' See \code{\link{sCorrect}} for more details.
#' @param ... arguments passed to the \code{summary} method of the object (when calling \code{summary})
#' or to the \code{sCorrect} method (when calling \code{summary2}). 
#' 
#' @seealso \code{\link{sCorrect}} for more detail about the small sample correction.
#'
#' @details \code{summary2} is the same as \code{summary}
#' except that it first computes the small sample correction (but does not store it).
#' So if \code{summary2} is to be called several times,
#' it is more efficient to pre-compute the quantities for the small sample correction
#' using \code{sCorrect} and then call \code{summary}.
#' 
#' @examples
#' m <- lvm(Y~X1+X2)
#' set.seed(10)
#' d <- sim(m, 2e1)
#'
#' ## Gold standard
#' summary(lm(Y~X1+X2, d))$coef
#' 
#' ## gls models
#' library(nlme)
#' e.gls <- gls(Y~X1+X2, data = d, method = "ML")
#' summary(e.gls)$tTable
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- FALSE ## no small sample correction
#' summary(e.gls)$tTable
#' 
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE ## small sample correction
#' summary(e.gls)$tTable
#' 
#' ## lvm models
#' e.lvm <- estimate(m, data = d)
#' summary(e.lvm)$coef
#' 
#' sCorrect(e.lvm) <- FALSE ## no small sample correction
#' summary(e.lvm)$coef
#' 
#' sCorrect(e.lvm) <- TRUE ## small sample correction
#' summary(e.lvm)$coef
#' 
#' @concept small sample inference

## * summary.gls2
#' @rdname summary
#' @method summary gls2
#' @export
summary.gls2 <- function(object, 
                         digit = max(3, getOption("digit")),
                         ...){
    
    ### ** perform Wald test
    name.param <- names(coef(object))
    n.param <- length(name.param)

    tTable.all <- compare2(object, par = name.param, as.lava = FALSE)
    tTable <- tTable.all[1:n.param,c("estimate","std","statistic","p-value","df")]
    dimnames(tTable) <- list(name.param,
                             c("Value","Std.Error","t-value","p-value","df")
                             )

    ### ** get summary
    class(object) <- setdiff(class(object),"gls2")
    object.summary <- summary(object, digits = digit, ...)
    
    ### ** update summary
    object.summary$tTable <- tTable

    ### ** export
    return(object.summary)
}

## * summary.lme2
#' @rdname summary
#' @method summary lme2
#' @export
summary.lme2 <- summary.gls2

## * summary.lvmfit2
#' @rdname summary
#' @method summary lvmfit2
#' @export
summary.lvmfit2 <- function(object, ...){

    ### ** perform Wald test
    param <- lava::pars(object)
    name.param <- names(param)
    n.param <- length(param)

    table.all <- compare2(object, par = name.param, as.lava = FALSE)
    table.coef <- table.all[1:n.param,c("estimate","std","statistic","p-value","df")]
    dimnames(table.coef) <- list(name.param,
                                 c("Estimate", "Std. Error", "t-value", "P-value", "df")
                                 )

    ### ** get summary
    class(object) <- setdiff(class(object),"lvmfit2")
    object.summary <- summary(object, ...)

    ## find digit
    vec.char <- setdiff(object.summary$coefmat[,"Estimate"],"")
    digit <- max(c(nchar(gsub(".","",vec.char,fixed = TRUE)))-1,1)

    ### ** update summary
    ### *** vcov
    object.summary$vcov <- attr(object$dVcov, "vcov.param")[name.param,name.param]    

    ### *** coef
    ## re-order table according to lava
    table.coef <- table.coef[rownames(object.summary$coef),,drop=FALSE]
    ## remove unappropriate p.values
    table.coef[is.na(object.summary$coef[,"P-value"]),"P-value"] <- NA
    object.summary$coef <- table.coef
    
    ### *** coefmat
    name.label0 <- trimws(rownames(CoefMat(object, labels = 0, level = 9)), which = "both")
    index.titleVariance <- which(name.label0=="Residual Variances:")
    if(length(index.titleVariance)>0){
        index.titleVariance <- (index.titleVariance+1):length(name.label0)
        name.label0[index.titleVariance] <- paste0(name.label0[index.titleVariance],lava.options()$symbols[2],name.label0[index.titleVariance])
    }

    table.coefmat <- object.summary$coefmat
    colnames(table.coefmat)[3:5] <- c("t-value","P-value","df")
    
    ## mimic lava:::CoefMat (called by lava:::summary.lvmfit)    
    e2add <- format(round(table.coef[,"Estimate"], max(1, digit - 1)), digits = digit - 1)
    e2add <- gsub(" NA","",e2add)
    sd2add <- format(round(table.coef[,"Std. Error"], max(1, digit - 1)), digits = digit - 1)
    sd2add <- gsub(" NA","",sd2add)
    df2add <- as.character(round(table.coef[,"df"],2))    
    df2add[is.na(df2add)] <- ""
    t2add <- format(round(table.coef[,"t-value"], max(1, digit - 1)), digits = digit - 1)
    t2add <- gsub(" NA","",t2add)

    p2add <- formatC(table.coef[,"P-value"], digits = digit - 1, format = "g",  preserve.width = "common", flag = "")
    p2add <- gsub(" NA","",p2add)
    p2add[table.coef[,"P-value"] < 1e-12] <- "  <1e-12"

    M2add <- cbind(e2add,sd2add,t2add,p2add,df2add)
    table.coefmat[match(rownames(table.coef), name.label0),] <- M2add

    table.coefmat[object.summary$coefma[,"P-value"]=="","P-value"] <- ""
    object.summary$coefmat <- table.coefmat

    ### ** Export
    return(object.summary)    
}

## * summary2
#' @rdname summary
#' @export
`summary2` <-
  function(object,...) UseMethod("summary2")

## * summary2.gls
#' @rdname summary
#' @export
summary2.gls <- function(object, bias.correct, ...){
    sCorrect(object, ...) <- bias.correct
    return(summary(object))
}

## * summary2.lme
#' @rdname summary
#' @export
summary2.lme <- summary2.gls

## * summary2.lvmfit
#' @rdname summary
#' @export
summary2.lvmfit <- function(object, bias.correct = TRUE, ...){
    sCorrect(object, ...) <- bias.correct
    return(summary(object))
}


##----------------------------------------------------------------------
### summary2.R ends here


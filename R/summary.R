### summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: jan 15 2018 (11:32) 
##           By: Brice Ozenne
##     Update #: 120
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
#' @param adjust.residuals Small sample correction: should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.
#' @param digit the number of digit to keep when diplaying the summary.
#' @param ... arguments passed to lower level methods.
#' 
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
#' dVcov2(e.gls, cluster = 1:NROW(d)) <- FALSE ## no small sample correction
#' summary(e.gls)$tTable
#' 
#' dVcov2(e.gls, cluster = 1:NROW(d)) <- TRUE ## small sample correction
#' summary(e.gls)$tTable
#' 
#' ## lvm models
#' e.lvm <- estimate(m, data = d)
#' summary(e.lvm)$coef
#' dVcov2(e.lvm) <- FALSE ## no small sample correction
#' summary(e.lvm)$coef
#' 
#' dVcov2(e.lvm) <- TRUE ## small sample correction
#' summary(e.lvm)$coef
#'

## * summary.gls2
#' @rdname summary
#' @method summary gls2
#' @export
summary.gls2 <- function(object, 
                         digit = max(3, getOption("digit")),
                         adjust.residuals = TRUE, ...){

    class(object) <- setdiff(class(object),"gls2")
    object.summary <- summary(object, digits = digit, ...)

    ## find digit
    
    
    ### ** update summary
    tTable <- lTest(object, Ftest = FALSE)[rownames(object.summary$tTable),c(1:3,5,4)]
    colnames(tTable) <- c("Value","Std.Error","t-value","p-value","df")

    object.summary$tTable <- tTable
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
summary.lvmfit2 <- function(object, adjust.residuals = FALSE, ...){

    class(object) <- setdiff(class(object),"lvmfit2")
    
    object.summary <- summary(object, ...)

    ## find digit
    vec.char <- setdiff(object.summary$coefmat[,"Estimate"],"")
    digit <- max(c(nchar(gsub(".","",vec.char,fixed = TRUE)))-1,1)

    ##
    param <- lava::pars(object)
    name.param <- names(param)
    name.allParam <- rownames(object.summary$coef)
    n.allParam <- length(name.allParam)
    data <- stats::model.frame(object)
    
    vcov.object <- attr(object$dVcov, "vcov.param")
        
    ### ** update summary
    ### *** vcov
    object.summary$vcov <- vcov.object[name.param,name.param]    

    ### *** coef
    table.coef <- lTest(object, Ftest = FALSE)[rownames(object.summary$coef),c(1:3,5,4)]
    colnames(table.coef) <- c("Estimate", "Std. Error", "t-value", "P-value", "df")
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

##----------------------------------------------------------------------
### summary2.R ends here

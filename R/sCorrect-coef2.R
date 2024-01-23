### sCorrect-coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:14) 
## Version: 
## Last-Updated: jan 23 2024 (10:23) 
##           By: Brice Ozenne
##     Update #: 321
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Model Coefficients With Small Sample Correction
#' @description Extract the coefficients from a latent variable model.
#' Similar to \code{lava::compare} but with small sample correction.
#' @name coef2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the model coefficients.
#' 
#' @return A numeric vector named with the names of the coefficients.
#' 
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#' 
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`coef2` <-
    function(object, as.lava, ...) UseMethod("coef2")


## * Examples
#' @rdname coef2
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' coef(e.lvm)
#' coef2(e.lvm)
#' coef2(e.lvm, as.lava = FALSE)

## * coef2.lvmfit
#' @export
coef2.lvmfit <- function(object, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(coef(estimate2(object, ssc = ssc, ...), as.lava = as.lava))

}

## * coef2.lvmfit2
#' @export
coef2.lvmfit2 <- function(object, as.lava = TRUE, ...){
    dots <- list(...)
    if(any(names(dots) %in% c("type","symbol","labels") == FALSE)){
        if(length(dots)>0){
            warning("Argument(s) \'",paste(setdiff(names(dots), c("type","symbol","labels")), collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
        }
    }
    if(length(dots)>1){ ## for the print function

        ## new values
        res <- model.tables(object, as.lava = TRUE)

        ## extract structure from lava
        object0 <- object
        class(object0) <- setdiff(class(object), "lvmfit2")
        out <- do.call(stats::coef, args = c(list(object0),dots)) ## this does not necessarily output the full parameter name
        ## Y1~~Y1 may be abreviated into Y1 which is confusion with Y1 the intercept
        dots$symbol <- NULL
        out.names <- intersect(rownames(do.call(stats::coef, args = c(list(object0),dots))), ## full name
                               rownames(res))
        index.out.names <- match(out.names, rownames(do.call(stats::coef, args = c(list(object0),dots))))

        ## rownames(res) <- as.character(object$sCorrect$skeleton$originalLink2param)
        out[index.out.names,"Estimate"] <- res[out.names,"estimate"]
        out[index.out.names,"Std. Error"] <- res[out.names,"se"]
        out[index.out.names,3] <- res[out.names,"statistic"] ## use 3 instead of Z value / Z-value
        if(object$sCorrect$df=="satterthwaite"){ 
            colnames(out)[3] <- "t value"
            if(colnames(out)[4]=="Pr(>|z|)"){ 
                colnames(out)[4] <- "Pr(>|t|)"
            }
        }
        out[index.out.names,4] <- res[out.names,"p.value"] ## use 4 instead of P-value / Pr(>|z|)

    }else{        
        out <- object$sCorrect$param[names(object$sCorrect$skeleton$originalLink2param)]
        if(as.lava==FALSE){
            names(out) <- as.character(object$sCorrect$skeleton$originalLink2param)
        }
    }
    return(out)
}

## * coef.lvmfit2
#' @export
coef.lvmfit2 <- coef2.lvmfit2

######################################################################
### sCorrect-coef2.R ends here

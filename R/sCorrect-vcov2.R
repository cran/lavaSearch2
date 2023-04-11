### sCorrect-vcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (13:55) 
## Version: 
## Last-Updated: jan 18 2022 (09:37) 
##           By: Brice Ozenne
##     Update #: 101
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov2 (documentation)
#' @title  Variance-Covariance With Small Sample Correction
#' @description  Extract the variance-covariance matrix from a latent variable model.
#' Similar to \code{stats::vcov} but with small sample correction.
#' @name vcov2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param robust [logical] should robust standard errors be used instead of the model based standard errors? Should be \code{TRUE} if argument cluster is not \code{NULL}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#' 
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the variance-covariance matrix.
#'
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#'
#' @return A matrix with as many rows and columns as the number of coefficients.
#' 
#' @examples
#' #### simulate data ####
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- Sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' #### linear models ####
#' e.lm <- lm(formula.lvm,data=d)
#'
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov0 <- vcov(e.lvm)
#' vcovSSC <- vcov2(e.lvm)
#' 
#' vcovSSC/vcov0
#' vcovSSC[1:4,1:4]/vcov(e.lm)
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`vcov2` <-
    function(object, robust, cluster, as.lava, ...) UseMethod("vcov2")

## * vcov2.lvmfit
#' @rdname vcov2
#' @export
vcov2.lvmfit <- function(object, robust = FALSE, cluster = NULL, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(vcov(estimate2(object, ssc = ssc, ...), robust = robust, cluster = cluster, as.lava = as.lava))

}

## * vcov2.lvmfit2
#' @rdname vcov2
#' @export
vcov2.lvmfit2 <- function(object, robust = FALSE, cluster = NULL, as.lava = TRUE, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }
    

    if(robust){
        out <- crossprod(iid(object, cluster = cluster, as.lava = as.lava))
    }else{
        if(!is.null(cluster)){
            warning("Argument \'cluster\' disregarded when argument \'robust\' is FALSE. \n")
        }
        out <- object$sCorrect$vcov.param
    }

    out <- out[names(object$sCorrect$skeleton$originalLink2param),
               names(object$sCorrect$skeleton$originalLink2param),
               drop=FALSE]
    if(as.lava==FALSE){
        dimnames(out) <- list(as.character(object$sCorrect$skeleton$originalLink2param),
                              as.character(object$sCorrect$skeleton$originalLink2param))
    }

    return(out)
}

## * vcov.lvmfit2
#' @rdname vcov2
#' @export
vcov.lvmfit2 <- vcov2.lvmfit2



######################################################################
### sCorrect-vcov2.R ends here

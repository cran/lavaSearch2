### mlf2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: jan 17 2018 (17:16) 
##           By: Brice Ozenne
##     Update #: 128
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mlf2
#' @title Simultaneous Inference for Multiple Models
#' @description Calculation of correlation between test statistics from multiple marginal models using the score decomposition. Similar to \code{multcomp::mmm} but with small sample corrections.
#'
#' @param ... A names argument list containing the fitted model.
#' See the documentation of \code{multcomp::mmm} for more details.
#'
#' @export
mlf2 <- function(...) {

     ret <- list(...)
     class(ret) <- "mlf2"
     ret

}

## * glht.mlf2
#' @title General Linear Hypothesis
#' @description Test general linear hypotheses and across latent variable models.
#' @name glht
#' 
#' @param model a list of latent variable models.
#' @param linfct a contrast matrix specifying the linear hypotheses to be tested.
#' @param adjust.residuals Small sample correction: should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.
#' @param robust should robust standard error be used? 
#' Otherwise rescale the influence function with the standard error obtained from the information matrix.
#' @param ... arguments passed to \code{glht}, \code{vcov}, and \code{lTest}.
#' 
#' @examples
#' library(multcomp)
#' 
#' ## Simulate data
#' mSim <- lvm(c(Y1,Y2,Y3)~ beta * eta, E ~ 1)
#' latent(mSim) <- "eta"
#' set.seed(10)
#' n <- 1e2
#'
#' df.data <- sim(mSim, n, latent = FALSE, p = c(beta = 1))
#'
#' ## Fit separate models
#' ls.lvm <- list(Y1 = estimate(lvm(Y1~E), data = df.data),
#'                Y2 = estimate(lvm(Y2~E), data = df.data),
#'                Y3 = estimate(lvm(Y3~E), data = df.data))
#' 
#' ## Create contrast matrix
#' C <- createContrast(ls.lvm, var.test = "E")
#'
#' lvm.glht <- glht2(ls.lvm, linfct = C)
#' summary(lvm.glht) ## adjusted
#' 
#' summary(lvm.glht, test = univariate()) ## not adjusted

#' @rdname glht
#' @method glht mlf2
#' @export 
glht.mlf2 <- function(model, linfct, ...) {

    ## ** check
    test.class <- c(inherits(model, c("mmm2")),
                    inherits(model, c("ls.lvmfit")))
    stopifnot(any(test.class)) ## necessary because only 2 methods defined to extract coef/vcov
    
    if (length(linfct) == 1) {
        linfct <- linfct[rep(1, length(model))]
        names(linfct) <- names(model)
    }
 
    ## ** get the name of each model
    n.model <- length(model)
    name.model <- names(model)    
    if(is.null(name.model)){
        stop("Argument ", sQuote("model")," must be named list. \n")
    }
        
    ## ** define contrast matrix
    if(is.list(linfct)){
        
        name.linfct <- names(linfct)
        if(is.null(name.linfct)){
            stop("Argument ", sQuote("linfct")," must be named list or a matrix. \n",
                 "Consider using the function createContrast when working with lvm objects. \n")
        }

        if (!identical(sort(name.linfct), sort(name.model))){
            stop("names of ", sQuote("model"), " and ", sQuote("linfct"),
                 " are not identical")
        }
    
        K <- lapply(name.model, function(i) {
            multcomp::glht(model[[i]], linfct = linfct[[i]], ...)$linfct
        })
    
        for (iK in 1:n.model) {
            rownames(K[[iK]]) <- paste(name.model[iK], rownames(K[[iK]]), sep = ": ")
            colnames(K[[iK]]) <- paste(name.model[iK], colnames(K[[iK]]), sep = ": ")
        }
        Cmatrix <- multcomp_.bdiag(K)
    }else{
        K <- lapply(name.model, function(x){ ## x <- name.model[2]
            iRownames <- grep(paste0(x,": "), rownames(linfct), value = FALSE, fixed = TRUE)
            iColnames <- grep(paste0(x,": "), colnames(linfct), value = FALSE, fixed = TRUE)
            linfct[iRownames, iColnames,drop=FALSE]            
        })
        Cmatrix <- linfct
    }
    names(K) <- name.model

    ## ** Get estimates
    ## do not add variance in the first run because cannot pass addtional arguments
    out <- multcomp_glht.matrix(model, linfct = Cmatrix, ...)

    ## ** Add variance
    ## call vcov.mmm2
    out$vcov <- stats::vcov(model, return.null = FALSE, ...)

    ## add degrees of freedom    
    df <- sapply(name.model, function(iName){ ## iName <- names(model)[1]
        if(identical(class(model[[iName]]),"lm") && "sigma2" %in% names(K[[iName]]) == FALSE){
            K[[iName]] <- cbind(K[[iName]], sigma2 = 0)
        }
        lTest(model[[iName]], C = K[[iName]], Ftest = FALSE, ...)$df
    })
    out$df <- as.double(round(stats::median(df)))
    return(out)
}

#' @rdname glht
#' @export
glht2 <- function (model, linfct, adjust.residuals = TRUE, robust = FALSE, ...){

    if(!is.matrix(linfct)){
        stop("Argument \'linfct\' must be a matrix \n")
    }

    if("lvmfit" %in% class(model)){
        out <- glht(model, linfct, ...)

        dVcov2(model, return.score = adjust.residuals) <- adjust.residuals
        res <- lTest(model, C = linfct, Ftest = TRUE)
        out$df <- round(res["global","df"])

        if(robust){
            out$vcov <- crossprod(attr(model$dVcov, "score") %*% attr(model$dVcov, "vcov.param"))
        }else{
            out$vcov <- attr(model$dVcov,"vcov.param")
        }

    }else{
        class(linfct) <- append("mlf2",class(linfct))
        class(model) <- "ls.lvmfit"
        out <- glht(model, linfct,
                    adjust.residuals = adjust.residuals, robust = robust, ...)
    }

    return(out)
}


##----------------------------------------------------------------------
### mlf2.R ends here

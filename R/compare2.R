### compare2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 30 2018 (14:33) 
## Version: 
## Last-Updated: jul 16 2018 (16:35) 
##           By: Brice Ozenne
##     Update #: 379
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - compare2
#' @title Test Linear Hypotheses with small sample correction
#' @description Test Linear Hypotheses using a multivariate Wald statistic.
#' Similar to \code{lava::compare} but with small sample correction.
#' @name compare2
#'
#' @param object an object that inherits from lm/gls/lme/lvmfit.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias? Argument passed to \code{sCorrect}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param par [vector of characters] expression defining the linear hypotheses to be tested.
#' See the examples section. 
#' @param contrast [matrix] a contrast matrix defining the left hand side of the linear hypotheses to be tested.
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param null,rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param as.lava [logical] should the output be similar to the one return by \code{lava::compare}?
#' @param F.test [logical] should a joint test be performed?
#' @param level [numeric 0-1] the confidence level of the confidence interval.
#' @param ...  [internal] only used by the generic method.
#'
#' @details The \code{par} argument or the arguments \code{contrast} and \code{null} (or equivalenty \code{rhs})
#' specify the set of linear hypotheses to be tested. They can be written:
#' \deqn{
#'   contrast * \theta = null
#' }
#' where \eqn{\theta} is the vector of the model coefficients. \cr
#' The \code{par} argument must contain expression(s) involving the model coefficients.
#' For example \code{"beta = 0"} or \code{c("-5*beta + alpha = 3","-alpha")} are valid expressions if alpha and beta belong to the set of model coefficients.
#' A contrast matrix and the right hand side will be generated inside the function. \cr
#' 
#' When directly specified, the contrast matrix must contain as many columns as there are coefficients in the model (mean and variance coefficients).
#' Each hypothesis correspond to a row in the contrast matrix. \cr
#'
#' The null vector should contain as many elements as there are row in the contrast matrix. \cr
#' 
#' Argument rhs and null are equivalent.
#' This redondance enable compatibility between \code{lava::compare}, \code{compare2}, \code{multcomp::glht}, and \code{glht2}.
#'
#' @seealso \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{sCorrect}} to pre-compute quantities for the small sample correction.
#' 
#' @return If \code{as.lava=TRUE} an object of class \code{htest}.
#' Otherwise a \code{data.frame} object.
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' mSim <- lvm(Y~0.1*X1+0.2*X2)
#' categorical(mSim, labels = c("a","b","c")) <- ~X1
#' transform(mSim, Id~Y) <- function(x){1:NROW(x)}
#' df.data <- lava::sim(mSim, 1e2)
#'
#' #### with lm ####
#' ## direct use of compare2
#' e.lm <- lm(Y~X1+X2, data = df.data)
#' anova(e.lm)
#' compare2(e.lm, par = c("X1b=0","X1c=0"))
#' 
#' ## or first compute the derivative of the information matrix
#' sCorrect(e.lm) <- TRUE
#' 
#' ## and define the contrast matrix
#' C <- createContrast(e.lm, par = c("X1b=0","X1c=0"), add.variance = TRUE)
#'
#' ## run compare2
#' compare2(e.lm, contrast = C$contrast, null = C$null)
#' compare2(e.lm, contrast = C$contrast, null = C$null, robust = TRUE)
#' 
#' #### with gls ####
#' library(nlme)
#' e.gls <- gls(Y~X1+X2, data = df.data, method = "ML")
#'
#' ## first compute the derivative of the information matrix
#' sCorrect(e.gls, cluster = 1:NROW(df.data)) <- TRUE
#' 
#' compare2(e.gls, par = c("5*X1b+2*X2 = 0","(Intercept) = 0"))
#' 
#' #### with lvm ####
#' m <- lvm(Y~X1+X2)
#' e.lvm <- estimate(m, df.data)
#' 
#' compare2(e.lvm, par = c("-Y","Y~X1b+Y~X1c"))
#' @concept small sample inference
#' @export
`compare2` <-
  function(object, ...) UseMethod("compare2")

## * compare2.lm
#' @rdname compare2
#' @export
compare2.lm <- function(object, df = TRUE, bias.correct = TRUE, ...){
    sCorrect(object, df = df) <- bias.correct
    return(.compare2(object, ...))
}

## * compare2.gls
#' @rdname compare2
#' @export
compare2.gls <- function(object, cluster = NULL, df = TRUE, bias.correct = TRUE, ...){
    sCorrect(object, df = df, cluster = cluster) <- bias.correct
    return(.compare2(object, ...))
}

## * compare2.lme
#' @rdname compare2
#' @export
compare2.lme <- compare2.lm

## * compare2.lvmfit
#' @rdname compare2
#' @export
compare2.lvmfit <- function(object, cluster = NULL, df = TRUE, bias.correct = TRUE, ...){
    sCorrect(object, df = df) <- bias.correct
    return(.compare2(object, cluster = cluster, ...))
}

## * compare2.lm2
#' @rdname compare2
#' @export
compare2.lm2 <- function(object, ...){
    return(.compare2(object, ...))
}

## * compare2.gls2
#' @rdname compare2
#' @export
compare2.gls2 <- function(object, ...){
    return(.compare2(object, ...))
}

## * compare2.lme2
#' @rdname compare2
#' @export
compare2.lme2 <- function(object, ...){
    return(.compare2(object, ...))
}

## * compare2.lvmfit2
#' @rdname compare2
#' @export
compare2.lvmfit2 <- function(object, ...){
    return(.compare2(object, ...))
}

## * .compare2
#' @rdname compare2
.compare2 <- function(object, par = NULL, contrast = NULL, null = NULL, rhs = NULL,
                      robust = FALSE, cluster = NULL, df = TRUE,
                      as.lava = TRUE, F.test = TRUE, level = 0.95){

    ## ** extract information
    if(!is.null(null) && !is.null(rhs)){
        stop("Arguments \'null\' and \'rhs\' should not be both specified \n")
    }

    if(robust && df > 1){ ## hidden option: df for robust standard error computed by (Pan, 2002)
        dVcov.param <- NA ## necessary to pass the test after
        score <- object$sCorrect$score
    }else if(df>0){        
        dVcov.param <- object$sCorrect$dVcov.param
    }else{
        dVcov.param <- NULL
    }
    
    param <- object$sCorrect$param
    if(robust){
        vcov.param <- crossprod(iid2(object, cluster = cluster))
    }else{
        vcov.param <- object$sCorrect$vcov.param
        attr(vcov.param, "warning") <- NULL
    }
    warn <- attr(object$sCorrect$vcov.param, "warning")
    keep.param <- dimnames(dVcov.param)[[3]]

    n.param <- length(param)
    name.param <- names(param)

    ### ** normalize linear hypotheses
    if(!is.null(par)){
        
        if(!is.null(contrast)){
            stop("Argument \'par\' and argument \'contrast\' should not simultaneously specified")
        }else if(!is.null(null)){
            stop("Argument \'par\' and argument \'null\' should not simultaneously specified")
        }else{
            res.C <- createContrast(par, name.param = name.param, add.rowname = TRUE)
            contrast <- res.C$contrast
            null <- res.C$null
        }
        
    }else{
        
        if(is.null(contrast)){
            stop("Argument \'contrast\' and argument \'par\' cannot be both NULL \n",
                 "Please specify the null hypotheses using one of the two arguments \n")
        }
        if(is.null(colnames(contrast))){
            stop("Argument \'contrast\' must have column names \n")
        }
        if(any(colnames(contrast) %in% name.param == FALSE)){
            txt <- setdiff(colnames(contrast), name.param)
            stop("Argument \'contrast\' has incorrect column names \n",
                 "invalid name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        if(any(name.param %in% colnames(contrast) == FALSE)){
            txt <- setdiff(name.param, colnames(contrast))
            stop("Argument \'contrast\' has incorrect column names \n",
                 "missing name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        if(NCOL(contrast) != n.param){
            stop("Argument \'contrast\' should be a matrix with ",n.param," columns \n")
        }
        ## reorder columns according to coefficients
        contrast <- contrast[,name.param,drop=FALSE]
        if(any(abs(svd(contrast)$d)<1e-10)){
            stop("Argument \'contrast\' is singular \n")
        }
        if(is.null(null)){
            null <- setNames(rep(0,NROW(contrast)),rownames(contrast))
        }else if(length(null)!=NROW(contrast)){
            stop("The length of argument \'null\' does not match the number of rows of argument \'contrast' \n")
        }
        if(is.null(rownames(contrast))){
            rownames(contrast) <- .contrast2name(contrast, null = null)
            null <- setNames(null, rownames(contrast))
        }
    }
    
    ### ** prepare export
    name.hypo <- rownames(contrast)
    n.hypo <- NROW(contrast)

    df.table <- as.data.frame(matrix(NA, nrow = n.hypo, ncol = 5,
                                     dimnames = list(name.hypo,
                                                     c("estimate","std","statistic","df","p-value"))
                                     ))

### ** Compute degrees of freedom

    if(is.null(dVcov.param)){
        df.Wald <- rep(Inf, n.hypo)
        df.F <- Inf
    }else{
        vcov.tempo <- object$sCorrect$vcov.param
        attr(vcov.tempo, "warning") <- NULL
        ## univariate
        if(robust && df > 1){ ## hidden option: df for robust standard error computed by (Pan, 2002)
            df.Wald  <- dfSigmaRobust(contrast = contrast,
                                      vcov = vcov.tempo,
                                      score = score)
        }else{
            df.Wald  <- dfSigma(contrast = contrast,
                                vcov = vcov.tempo,
                                dVcov = dVcov.param,
                                keep.param = keep.param)
        }

        ## multivariate
        svd.tempo <- eigen(solve(contrast %*% vcov.tempo %*% t(contrast)))
        D.svd <- diag(svd.tempo$values, nrow = n.hypo, ncol = n.hypo)
        P.svd <- svd.tempo$vectors
     
        C.anova <- sqrt(D.svd) %*% t(P.svd) %*% contrast
        ## Fstat - crossprod(C.anova %*% p)/n.hypo
        if(robust && df > 1){ ## hidden option: df for robust standard error computed by (Pan, 2002)
            nu_m  <- dfSigmaRobust(contrast = C.anova,
                                   vcov = vcov.tempo,
                                   score = score)
        }else{
            nu_m <- dfSigma(contrast = C.anova,
                            vcov = vcov.tempo,
                            dVcov = dVcov.param,
                            keep.param = keep.param) ## degree of freedom of the independent t statistics
        }    
        EQ <- sum(nu_m/(nu_m-2))
        df.F <- 2*EQ / (EQ - n.hypo)

    }

### ** Wald test
    ## statistic
    C.p <- (contrast %*% param) - null
    C.vcov.C <- contrast %*% vcov.param %*% t(contrast)
    sd.C.p <- sqrt(diag(C.vcov.C))
    stat.Wald <- C.p/sd.C.p
    
    ## store
    df.table$estimate <- as.numeric(C.p)
    df.table$std <- as.numeric(sd.C.p)
    df.table$statistic <- as.numeric(stat.Wald)
    df.table$df <- as.numeric(df.Wald)
    df.table$`p-value` <- as.numeric(2*(1-stats::pt(abs(df.table$statistic), df = df.table$df)))
    
### ** multivariate F test
    df.table <- rbind(df.table, global = rep(NA,5))
    error <- NULL
     
    if(F.test){
        ## statistic
        stat.F <- try(t(C.p) %*% solve(C.vcov.C)%*% (C.p) / n.hypo, silent = TRUE)
     
        ## store
        if(!inherits(stat.F,"try-error")){
            df.table["global", "statistic"] <- as.numeric(stat.F)
            df.table["global", "df"] <- df.F
            df.table["global", "p-value"] <- 1 - stats::pf(df.table["global", "statistic"],
                                                           df1 = n.hypo,
                                                           df2 = df.table["global", "df"])
        }else{
            error <- df.table
        }
    }
    
    ## ** export
    if(as.lava == TRUE){
        level.inf <- (1-level)/2
        level.sup <- 1-level.inf

        level.inf.label <- paste0(100*level.inf,"%")
        level.sup.label <- paste0(100*level.sup,"%")

        df.estimate <- matrix(NA, nrow = n.hypo, ncol = 5,
                              dimnames = list(name.hypo,c("Estimate", "Std.Err", "df", level.inf.label, level.sup.label)))
        df.estimate[,"Estimate"] <- df.table[name.hypo,"estimate"]
        df.estimate[,"Std.Err"] <- df.table[name.hypo,"std"]
        df.estimate[,"df"] <- df.table[name.hypo,"df"]
        df.estimate[,level.inf.label] <- df.table[name.hypo,"estimate"] + stats::qt(level.inf, df = df.table[name.hypo,"df"]) * df.table[name.hypo,"std"]
        df.estimate[,level.sup.label] <- df.table[name.hypo,"estimate"] + stats::qt(level.sup, df = df.table[name.hypo,"df"]) * df.table[name.hypo,"std"]

        out <- list(statistic = setNames(df.table["global","statistic"],"F-statistic"),
                    parameter = setNames(round(df.table["global","df"],2), paste0("df1 = ",n.hypo,", df2")), ## NOTE: cannot not be change to coefficients because of lava
                    p.value = df.table["global","p-value"],
                    method = c("- Wald test -", "", "Null Hypothesis:", name.hypo),
                    estimate = df.estimate,
                    vcov = C.vcov.C,
                    coef = C.p[,1],
                    null = null,
                    cnames = name.hypo                    
                    )
        if(robust){
            colnames(out$estimate)[2] <- "robust SE"
        }        
        attr(out, "B") <- contrast
        class(out) <- "htest"
    }else{
        out <- df.table
        attr(out, "warning") <- warn
        attr(out, "contrast") <- contrast
    }

    attr(out,"error") <- error
    return(out)
}

## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param vcov [numeric matrix] the variance-covariance matrix of the parameters.
##' @param dVcov [numeric array] the first derivative of the variance-covariance matrix of the parameters.
##' @param keep.param [character vector] the name of the parameters with non-zero first derivative of their variance parameter.
##' 
dfSigma <- function(contrast, vcov, dVcov, keep.param){
    C.vcov.C <- rowSums(contrast %*% vcov * contrast) ## variance matrix of the linear combination
    
    C.dVcov.C <- sapply(keep.param, function(x){
        rowSums(contrast %*% dVcov[,,x] * contrast)
    })
    numerator <- 2 *(C.vcov.C)^2
    denom <- rowSums(C.dVcov.C %*% vcov[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
    df <- numerator/denom
    return(df)
}

## * dfSigmaRobust
##' @title Degree of Freedom for the Robust Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the robust-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param vcov [numeric matrix] the variance-covariance matrix of the parameters.
##' @param score [numeric matrix] the individual score for each parameter.
##'
##' @details When contrast is the identity matrix, this function compute the moments of the sandwich estimator
##' and the degrees of freedom of the approximate t-test as described in (Pan, 2002) section 2 and 3.1.
##'
##' @references
##' Wei Pan and Melanie M. Wall, Small-sample adjustments in using the sandwich variance estiamtor in generalized estimating equations. Statistics in medicine (2002) 21:1429-1441.
##' 
dfSigmaRobust <- function(contrast, vcov, score){
    
    ## ** prepare
    n <- NROW(score)
    p <- NCOL(contrast)
        
    ## apply contrasts
    score.S <- score %*% contrast
    vcov.S <- vcov %*% contrast
    
    M.vecPi <- t(apply(score.S, 1, function(iRow){
        as.vector(tcrossprod(iRow))
    }))
    ## dim(M.vecPi) ## (n,p)

    #### ** check
    ## vcov %*% t(score) %*% score %*% vcov
    ## (vcov %x% vcov) %*% colSums(M.vecPi)

    ## ** compute moments of P
    Q <- colMeans(M.vecPi)
    T <- stats::var(M.vecPi) ## missing 1/n factor compared to (Pan, 2002)
    
    ## vec.Pcenter <- sweep(vec.P, FUN = "-", MARGIN = 2, STATS = Q)
    ## crossprod(vec.Pcenter)/(99)

    ## ** compute moments of Vs = (Vm x Vm) P
    E_Vs <- n *(vcov.S %x% vcov.S) %*% Q
    Sigma_Vs <- n * (vcov.S %x% vcov.S) %*% T %*% (vcov.S %x% vcov.S) ## missing n factor compared to (Pan, 2002) so it cancels out with the previous missing factor

    ## ** compute parameters of the chi-square distribution
    sigma <- as.double(E_Vs)
    tau <- diag(Sigma_Vs)

    df <- matrix((2*sigma^2)/tau, nrow = p, ncol = p)
    return(diag(df))
}


##----------------------------------------------------------------------
### compare2.R ends here

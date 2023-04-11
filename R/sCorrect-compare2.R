### compare2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 30 2018 (14:33) 
## Version: 
## Last-Updated: Apr 11 2023 (22:31) 
##           By: Brice Ozenne
##     Update #: 903
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - compare2
#' @title Test Linear Hypotheses With Small Sample Correction
#' @description Test Linear Hypotheses using Wald statistics in a latent variable model.
#' Similar to \code{lava::compare} but with small sample correction.
#' @name compare2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param linfct [matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] should the output be similar to the one return by \code{lava::compare}?
#' @param F.test [logical] should a joint test be performed?
#' @param conf.level [numeric 0-1] level of the confidence intervals.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/code{FALSE}/code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details The \code{linfct} argument and \code{rhs} specify the set of linear hypotheses to be tested. They can be written:
#' \deqn{
#'   linfct * \theta = rhs
#' }
#' where \eqn{\theta} is the vector of the model coefficients. \cr
#' The \code{par} argument must contain expression(s) involving the model coefficients.
#' For example \code{"beta = 0"} or \code{c("-5*beta + alpha = 3","-alpha")} are valid expressions if alpha and beta belong to the set of model coefficients.
#' A contrast matrix and the right hand side will be generated inside the function. \cr
#' 
#' When directly specified, the contrast matrix must contain as many columns as there are coefficients in the model (mean and variance coefficients).
#' Each hypothesis correspond to a row in the contrast matrix. \cr
#'
#' The rhs vector should contain as many elements as there are row in the contrast matrix. \cr
#' 
#' @seealso \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#' 
#' @return If \code{as.lava=TRUE} an object of class \code{htest}.
#' Otherwise a \code{data.frame} object.

## * example - compare2
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' mSim <- lvm(Y~0.1*X1+0.2*X2)
#' categorical(mSim, labels = c("a","b","c")) <- ~X1
#' transform(mSim, Id~Y) <- function(x){1:NROW(x)}
#' df.data <- lava::sim(mSim, 1e2)
#'
#' #### with lvm ####
#' m <- lvm(Y~X1+X2)
#' e.lvm <- estimate(m, df.data)
#' 
#' compare2(e.lvm, linfct = c("Y~X1b","Y~X1c","Y~X2"))
#' compare2(e.lvm, linfct = c("Y~X1b","Y~X1c","Y~X2"), robust = TRUE)
#' 
#' @concept inference
#' @keywords smallSampleCorrection
#' @export
`compare2` <-
    function(object, linfct, rhs,
             robust, cluster,
             as.lava, F.test,
             conf.level, ...) UseMethod("compare2")

## * compare2.lvmfit
#' @rdname compare2
#' @export
compare2.lvmfit <- function(object, linfct = NULL, rhs = NULL,
                            robust = FALSE, cluster = NULL,
                            as.lava = TRUE, F.test = TRUE,
                            conf.level = 0.95,
                            ssc = lava.options()$ssc, df = lava.options()$df, ...){

    return(compare(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...),
                    linfct = linfct, rhs = rhs, robust = robust, cluster = cluster, as.lava = as.lava, F.test = F.test, conf.level = conf.level)
           )

}

## * compare2.lvmfit2
#' @rdname compare2
#' @export
compare2.lvmfit2 <- function(object, linfct = NULL, rhs = NULL,
                              robust = FALSE, cluster = NULL,
                              as.lava = TRUE, F.test = TRUE,
                              conf.level = 0.95, ...){
   
    dots <- list(...)
    if(any(names(dots)=="par")){
        stop("Argument \'par\' is no longer used as it has been replaced by \'linfct\'. \n")
    }
    if(length(dots[names(dots) %in% "sep" == FALSE])>0){
        warning("Argument(s) \'",paste(setdiff(names(dots),"sep"),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }
    if(is.null(linfct)){ ## necessary for lava::gof to wo
        object0 <- object
        class(object0) <- setdiff(class(object0),"lvmfit2")
        return(lava::compare(object0))
    }
    if(!is.logical(robust)){ 
        stop("Argument \'robust\' should be TRUE or FALSE \n")
    }
    if(robust==FALSE && !is.null(cluster)){
        stop("Argument \'cluster\' must be NULL when argument \'robust\' is FALSE \n")
    }

    ## ** extract information
    df <- object$sCorrect$df
    
    ## 0-order: param
    param <- coef(object, as.lava = FALSE)
    n.param <- length(param)
    name.param <- names(param)

    ## 1-order: score
    if(robust){
        score <- score(object, cluster = cluster, as.lava = FALSE, indiv = TRUE)
    }else{
        score <- NULL
    }
    
    ## 2-order: variance covariance
    vcov.param <- vcov(object, as.lava = FALSE)
    warn <- attr(vcov.param, "warning")
    attr(vcov.param, "warning") <- NULL
    if(robust){
        rvcov.param <- vcov(object, robust = TRUE, cluster = cluster, as.lava = FALSE)
    }

    ## 3-order: derivative of the variance covariance matrices
    if(df == "satterthwaite"){
        dVcov.param <- object$sCorrect$dVcov.param[names(object$sCorrect$skeleton$originalLink2param),
                                                   names(object$sCorrect$skeleton$originalLink2param),
                                                   names(object$sCorrect$skeleton$originalLink2param),
                                                   drop=FALSE]
        dimnames(dVcov.param) <- list(as.character(object$sCorrect$skeleton$originalLink2param),
                                      as.character(object$sCorrect$skeleton$originalLink2param),
                                      as.character(object$sCorrect$skeleton$originalLink2param))
        keep.param <- dimnames(dVcov.param)[[3]]

        if(robust && (lava.options()$df.robust != 1)){
            if(!is.null(cluster) || is.null(object$sCorrect$dRvcov.param)){
                ## update derivative according to cluster
                hessian <- hessian2(object, cluster = cluster, as.lava = FALSE)
                dRvcov.param <- .dRvcov.param(score = score,
                                              hessian = hessian,
                                              vcov.param = vcov.param,
                                              dVcov.param = dVcov.param,
                                              n.param = n.param,
                                              name.param = name.param)
                                              
            }else{
                dRvcov.param <- object$sCorrect$dRvcov.param[names(object$sCorrect$skeleton$originalLink2param),
                                                             names(object$sCorrect$skeleton$originalLink2param),
                                                             names(object$sCorrect$skeleton$originalLink2param),
                                                             drop=FALSE]
                dimnames(dRvcov.param) <- list(as.character(object$sCorrect$skeleton$originalLink2param),
                                               as.character(object$sCorrect$skeleton$originalLink2param),
                                               as.character(object$sCorrect$skeleton$originalLink2param))
            }
        }
    }

    ## ** normalize linear hypotheses
    if(!is.matrix(linfct)){
        res.C <- createContrast(object, linfct = linfct, rowname.rhs = FALSE, ...)
        if(any(colnames(res.C$contrast)!=name.param) && all(colnames(res.C$contrast) == names(object$sCorrect$skeleton$originalLink2param))){
            colnames(res.C$contrast) <- as.character(object$sCorrect$skeleton$originalLink2param)
        }
        linfct <- res.C$contrast
        if(is.null(rhs)){
            rhs <- res.C$null
        }else{
            if(length(rhs)!=length(res.C$null)){
                stop("Incorrect argument \'rhs\' \n",
                     "Must have length ",length(res.C$null),"\n")
            }
            rhs <- stats::setNames(rhs, names(res.C$null))
        }
        name.hypoShort <- rownames(linfct)
        name.hypo <- paste0(name.hypoShort," = ",rhs)
    }else{
        if(is.null(colnames(linfct))){
            stop("Argument \'linfct\' must have column names \n")
        }
        if(NCOL(linfct) != n.param){
            stop("Argument \'linfct\' should be a matrix with ",n.param," columns \n")
        }
        if(any(colnames(linfct) %in% name.param == FALSE)){
            txt <- setdiff(colnames(linfct), name.param)
            stop("Argument \'linfct\' has incorrect column names \n",
                 "invalid name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        if(any(name.param %in% colnames(linfct) == FALSE)){
            txt <- setdiff(name.param, colnames(linfct))
            stop("Argument \'linfct\' has incorrect column names \n",
                 "missing name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        ## reorder columns according to coefficients
        linfct <- linfct[,name.param,drop=FALSE]
        if(F.test && any(abs(svd(linfct)$d)<1e-10)){
            stop("Argument \'linfct\' is singular \n")
        }
        if(is.null(rhs)){
            rhs <- stats::setNames(rep(0,NROW(linfct)),rownames(linfct))
        }else if(length(rhs)!=NROW(linfct)){
            stop("The length of argument \'rhs\' must match the number of rows of argument \'linfct' \n")
        }
        if(is.null(rownames(linfct))){
            rownames(linfct) <- .contrast2name(linfct, null = rhs)
            rhs <- stats::setNames(rhs, rownames(linfct))
        }
        name.hypo <- rownames(linfct)
        name.hypoShort <- sapply(strsplit(name.hypo, split = " = ", fixed = TRUE),"[[",1)
    }

    n.hypo <- length(name.hypo)
    linfct <- linfct[,names(param),drop=FALSE] ## column in contrast may not be in the same order as param

    ## ** Univariate Wald test
    ## coefficient (used for F.test and lava export)
    C.p <- linfct %*% param
    C.p.rhs <- C.p - rhs

    ## variance (used for F.test and lava export)
    if(robust){
        C.vcov.C <- linfct %*% rvcov.param %*% t(linfct)
    }else{
        C.vcov.C <- linfct %*% vcov.param %*% t(linfct)
    }

    ## df
    if(df == "satterthwaite"){
        df.Wald  <- dfSigma(contrast = linfct,
                            score = score,
                            vcov = vcov.param,
                            rvcov = rvcov.param,
                            dVcov = dVcov.param,
                            dRvcov = dRvcov.param,
                            keep.param = keep.param,                            
                            type = if(robust){lava.options()$df.robust}else{1})

        ##
        ## 2 * vcov.param["Y","Y"]^2 / (vcov.param["Y~~Y","Y~~Y"]*dVcov.param["Y","Y","Y~~Y"]^2)
        ## 
    }else{
        df.Wald <- rep(Inf, n.hypo)
    }

    ## ** Multivariate Wald test
    error <- NULL
    if(F.test){
        iC.vcov.C <- try(solve(C.vcov.C), silent = TRUE)
        if(!inherits(iC.vcov.C,"try-error")){
            stat.F <- t(C.p.rhs) %*% iC.vcov.C %*% (C.p.rhs) / n.hypo

            ## df (independent t statistics)
            if(df == "satterthwaite"){
                svd.tempo <- eigen(iC.vcov.C)
                D.svd <- diag(svd.tempo$values, nrow = n.hypo, ncol = n.hypo)
                P.svd <- svd.tempo$vectors
     
                C.anova <- sqrt(D.svd) %*% t(P.svd) %*% linfct

                nu_m  <- dfSigma(contrast = C.anova,
                                 score = score,
                                 vcov = vcov.param,
                                 rvcov = rvcov.param,
                                 dVcov = dVcov.param,
                                 dRvcov = dRvcov.param,
                                 keep.param = keep.param,                            
                                 type = if(robust){lava.options()$df.robust}else{1})
                EQ <- sum(nu_m/(nu_m-2))
                df.F <- 2*EQ / (EQ - n.hypo)
            }else{
                df.F <- Inf
            }
            ## store
            F.res <- c("statistic" = as.numeric(stat.F),
                       "df" = df.F,
                       "p.value" = 1 - stats::pf(stat.F,
                                                 df1 = n.hypo,
                                                 df2 = df.F)
                       )
        }else{
            warning("Unable to invert the variance-covariance matrix after application of the contrasts \n")
            error <- iC.vcov.C
        }
    }

    ## ** export
    if(as.lava == TRUE){
        level.inf <- (1-conf.level)/2
        level.sup <- 1-level.inf

        level.inf.label <- paste0(100*level.inf,"%")
        level.sup.label <- paste0(100*level.sup,"%")

        df.estimate <- matrix(NA, nrow = n.hypo, ncol = 5,
                              dimnames = list(name.hypoShort,c("Estimate", "Std.Err", "df", level.inf.label, level.sup.label)))
        df.estimate[,"Estimate"] <- C.p
        df.estimate[,"Std.Err"] <- sqrt(diag(C.vcov.C))
        df.estimate[,"df"] <- df.Wald
        df.estimate[,level.inf.label] <- df.estimate[,"Estimate"] + stats::qt(level.inf, df = df.estimate[,"df"]) * df.estimate[,"Std.Err"]
        df.estimate[,level.sup.label] <- df.estimate[,"Estimate"] + stats::qt(level.sup, df = df.estimate[,"df"]) * df.estimate[,"Std.Err"]

        dimnames(C.vcov.C) <- list(name.hypoShort,name.hypoShort)
        out <- list(statistic = stats::setNames(F.res["statistic"],"F-statistic"),
                    parameter = stats::setNames(round(F.res["df"],2), paste0("df1 = ",n.hypo,", df2")), ## NOTE: cannot not be change to coefficients because of lava
                    p.value = F.res["p.value"],
                    method = c("- Wald test -", "", "Null Hypothesis:", name.hypo),
                    estimate = df.estimate,
                    vcov = C.vcov.C,
                    coef = stats::setNames(C.p[,1], name.hypoShort),
                    null = stats::setNames(rhs, name.hypoShort),
                    cnames = name.hypo                    
                    )
        if(robust){
            colnames(out$estimate)[2] <- "robust SE"
        }
        rownames(linfct) <- name.hypo
        attr(out, "B") <- linfct
        class(out) <- "htest"
    }else{
        if(length(unique(df.Wald))==1){
            df.Wald <- df.Wald[1]
        }
        out <- list(model = object,
                    linfct = linfct,
                    rhs = unname(rhs),
                    coef = param,
                    vcov = if(robust){rvcov.param}else{vcov.param},
                    df = df.Wald,
                    alternative = "two.sided",
                    type = NULL,
                    robust = robust,
                    ssc = object$sCorrect$ssc$type,
                    global = if(F.test){F.res}else{NULL})
        class(out) <- c("glht2","glht")
    }
    attr(out,"warning") <- warn
    attr(out,"error") <- error
    return(out)

}

## * compare.lvmfit2
#' @rdname compare2
#' @export
compare.lvmfit2 <- compare2.lvmfit2

## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param score [numeric matrix] the individual score for each parameter.
##' @param vcov [numeric matrix] the model-based variance-covariance matrix of the parameters.
##' @param rvcov [numeric matrix] the robust variance-covariance matrix of the parameters.
##' @param dVcov [numeric array] the first derivative of the model-based variance-covariance matrix of the parameters.
##' @param dRvcov [numeric array] the first derivative of the robust variance-covariance matrix of the parameters.
##' @param keep.param [character vector] the name of the parameters with non-zero first derivative of their variance parameter.
##' @param type [integer] 1 corresponds to the Satterthwaite approximation of the the degrees of freedom applied to the model-based variance,
##' 2 to the Satterthwaite approximation of the the degrees of freedom applied to the robust variance,
##' 3 to the approximation described in (Pan, 2002) section 2 and 3.1.
##'
##' @references
##' Wei Pan and Melanie M. Wall, Small-sample adjustments in using the sandwich variance estiamtor in generalized estimating equations. Statistics in medicine (2002) 21:1429-1441.
##' 
dfSigma <- function(contrast, score, vcov, rvcov, dVcov, dRvcov, keep.param, type){
    if(type==1){
        C.vcov.C <- rowSums(contrast %*% vcov * contrast) ## variance matrix of the linear combination
        C.dVcov.C <- sapply(keep.param, function(x){
            rowSums(contrast %*% dVcov[,,x] * contrast)
        })
        numerator <- 2 *(C.vcov.C)^2
        denom <- rowSums(C.dVcov.C %*% vcov[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
        df <- numerator/denom
    }else if(type==2){
        C.rvcov.C <- rowSums(contrast %*% rvcov * contrast) ## variance matrix of the linear combination
        C.dRvcov.C <- sapply(keep.param, function(x){
            rowSums(contrast %*% dRvcov[,,x] * contrast)
        })
        numerator <- 2 *(C.rvcov.C)^2
        denom <- rowSums(C.dRvcov.C %*% rvcov[keep.param,keep.param,drop=FALSE] * C.dRvcov.C)
        df <- numerator/denom
    }else if(type==3){
        vcov.S <- contrast %*% vcov
        index.var <- diag(matrix(1:NROW(contrast)^2,NROW(contrast),NROW(contrast)))
        
        K <- NROW(score)
        ls.Pi <- lapply(1:K, function(iC){as.double(tcrossprod(score[iC,]))})
        M.Pi <- do.call(rbind,ls.Pi)
        M.Pi_center <- sweep(M.Pi, MARGIN = 2, STATS = colMeans(M.Pi), FUN = "-")
        ## M.Pi_center - M.Pi
        T <- t(M.Pi_center) %*% M.Pi_center / (K*(K-1))
        ## range(var(M.Pi)/K - T)
        eq.3 <- K^2 * (vcov.S %x% vcov.S) %*% T %*%  (vcov.S %x% vcov.S)

        Vs <- (vcov.S %x% vcov.S) %*% Reduce("+",ls.Pi)
        ## range(Vs - as.double(rvcov))
        ## range(Vs[index.var] - diag(rvcov))
        df <- 2*Vs[index.var]^2/sapply(index.var, function(iIndex){eq.3[iIndex,iIndex]})
        
    }
    
    return(stats::setNames(df, rownames(contrast)))
}


##----------------------------------------------------------------------
### compare2.R ends here

### lTest.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: jan 19 2018 (11:12) 
##           By: Brice Ozenne
##     Update #: 637
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - lTest
#' @title  Compute the Degree of Freedom of the Variance Parameters
#' @description Compute the degree of freedom of the variance parameters.
#' @name lTest
#'
#' @param object a lvm object.
#' @param C [optional] a contrast matrix.
#' @param adjust.residuals should a small-sample correction be used when computing the variance of the parameters and the degree of freedoms.
#' @param Ftest should a join test be computed.
#' @param ... arguments to be passed to \code{\link{dVcov2}}.
#'
#' @details In the case of a \code{lm} object, the contrast matrix need not to contain
#' a column for the variance parameter when the columns are named. When so, a column containing 0 is added to the contrast matrix.
#' 
#' @examples
#' set.seed(10)
#' mSim <- lvm(Y~0.1*X1+0.2*X2)
#' categorical(mSim, labels = c("a","b","c")) <- ~X1
#' transform(mSim, Id~Y) <- function(x){1:NROW(x)}
#' df.data <- lava::sim(mSim, 1e2)
#'
#' ## gold standard
#' e.lm <- lm(Y~X1+X2, data = df.data)
#' anova(e.lm)
#' 
#' lTest(e.lm)
#' 
#' ## gls model
#' library(nlme)
#' e.gls <- gls(Y~X1+X2, data = df.data, method = "ML")
#' e.gls$dVcov <- dVcov2(e.gls, data = df.data, cluster = df.data$Id)
#' 
#' C <- rbind(c(0,1,0,0,0),c(0,0,1,0,0))
#' colnames(C) <- names(attr(e.gls$dVcov,"param"))
#' lTest(e.gls, data = df.data, C = C)
#'
#' C <- rbind(c(0,0,0,1,0))
#' colnames(C) <- names(attr(e.gls$dVcov,"param"))
#' lTest(e.gls, data = df.data, C = C)
#' 
#' ## latent variable model
#' m <- lvm(Y~X1+X2)
#' e.lvm <- estimate(m, df.data)
#' e.lvm$dVcov <- dVcov2(e.lvm)
#' 
#' C <- rbind(c(0,0,1,0,0),c(0,0,0,1,0))
#' colnames(C) <- names(coef(e.lvm))
#' lTest(e.lvm, C = C)
#'
#' C <- rbind(c(0,1,0,0,0))
#' colnames(C) <- names(coef(e.lvm))
#' lTest(e.lvm, C = C)
#' 
#' @export
`lTest` <-
  function(object, ...) UseMethod("lTest")

## * lTest.lm
#' @rdname lTest
#' @export
lTest.lm <- function(object, C = NULL, adjust.residuals = TRUE,
                     Ftest = TRUE, ...){

    ## ** Extract information
    if(is.null(object$dVcov)){
        dVcov.dtheta  <- dVcov2(object, adjust.residuals = adjust.residuals, ...)
    }else{
        dVcov.dtheta <- object$dVcov
    }
    p <- attr(dVcov.dtheta, "param")

    vcov.param <- attr(dVcov.dtheta, "vcov.param")
    warn <- attr(vcov.param, "warning")
    attr(dVcov.dtheta, "vcov.param") <- NULL
    keep.param <- dimnames(dVcov.dtheta)[[3]]

    n.param <- length(p)
    name.param <- names(p)

    ### ** normalize C
    if(is.null(C)){
        C <- diag(1, nrow = n.param, ncol = n.param)
        dimnames(C) <- list(name.param, name.param)
    }else{
        if(NCOL(C) == (n.param-1) && all(names(stats::coef(object)) %in% colnames(C)) ){
            C <- cbind(C, sigma = 0)
            }else if(NCOL(C) != n.param){
            stop("Argument \'C\' should be a matrix with ",n.param," columns \n")
        }
        ## if(is.null(colnames(C)) || any(colnames(C) != name.param)){
        ##     stop("Argument \'C\' has incorrect column names \n")
        ## }
        if(any(abs(svd(C)$d)<1e-10)){
            stop("Argument \'C\' is singular \n")
        }
        if(is.null(colnames(C))){
            colnames(C) <- name.param            
        }
        if(is.null(rownames(C))){
            rownames(C) <- .contrast2name(C)
        }
    }
    q <- NROW(C)
       
    ### ** Compute degrees of freedom
    df.table <- as.data.frame(matrix(NA, nrow = q, ncol = 5))
    colnames(df.table) <- c("estimate","std","statistic","df","p-value")
    rownames(df.table) <- rownames(C)

    calcDF <- function(M.C){ # M.C <- C
        C.vcov.C <- rowSums(M.C %*% vcov.param * M.C)
    
        C.dVcov.C <- sapply(keep.param, function(x){
            rowSums(M.C %*% dVcov.dtheta[,,x] * M.C)
        })
        
        numerator <- 2 *(C.vcov.C)^2
        denom <- rowSums(C.dVcov.C %*% vcov.param[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
        df <- numerator/denom
        return(df)
    }

    ### *** Wald test
    ## statistic
    C.p <- C %*% p
    C.vcov.C <- C %*% vcov.param %*% t(C)
    sd.C.p <- sqrt(diag(C.vcov.C))
    stat.Wald <- C.p/sd.C.p
    
    ## df
    df.Wald  <- calcDF(C)
    ## if(adjust.residuals==FALSE){ ## for univariate linear models
    ##     df[name.coef] <- n
    ##     df["sigma2"] <- n/4
    ## }else{
    ##     df[name.coef] <- n^2/(n+p)
    ##     df["sigma2"] <- n^2/(4*(n+p))
    ## }
    
    ## store
    df.table[rownames(C), "estimate"] <- C.p
    df.table[rownames(C), "std"] <- sd.C.p
    df.table[rownames(C), "statistic"] <- stat.Wald
    df.table[rownames(C), "df"] <- df.Wald
    df.table[rownames(C), "p-value"] <- 2*(1-stats::pt(abs(df.table[rownames(C), "statistic"]),
                                                df = df.table[rownames(C), "df"]))
    
    ### *** F test
    if(Ftest){
        i.C.vcov.C <- solve(C.vcov.C)
        stat.F <- t(C.p) %*% i.C.vcov.C %*% (C.p) / q

        ## df
        svd.tempo <- eigen(i.C.vcov.C)
        D.svd <- diag(svd.tempo$values, nrow = q, ncol = q)
        P.svd <- svd.tempo$vectors
     
        C.anova <- sqrt(D.svd) %*% t(P.svd) %*% C
        ## Fstat - crossprod(C.anova %*% p)/q
        nu_m <- calcDF(C.anova) ## degree of freedom of the independent t statistics
    
        EQ <- sum(nu_m/(nu_m-2))
        df.F <- 2*EQ / (EQ - q)

        ## store
        df.table <- rbind(df.table, global = c(NA,NA,NA,NA,NA))
        df.table["global", "statistic"] <- as.numeric(stat.F)
        df.table["global", "df"] <- df.F
        df.table["global", "p-value"] <- 1 - stats::pf(df.table["global", "statistic"],
                                                       df1 = q,
                                                       df2 = df.table["global", "df"])
    }
    
    ## ** export
    attr(df.table, "warning") <- warn
    return(df.table)
}

## * lTest.gls
#' @rdname lTest
#' @export
lTest.gls <- lTest.lm

## * lTest.lme
#' @rdname lTest
#' @export
lTest.lme <- lTest.lm

## * lTest.lvmfit
#' @rdname lTest
#' @export
lTest.lvmfit <- lTest.lm

## * .contrast2name
#' @title Create Rownames for a Contrast Matrix
#' @description Create rownames for a contrast matrix.
#' using the coefficients and the names of the parameters
#' @name contrast2name
#'
#' @param C a constrast matrix.
#' 
.contrast2name <- function(C){
    col <- coef <- coefname <- rowname <- nrow <- NULL

    df.index <- as.data.frame(which(C != 0, arr.ind = TRUE))
    df.index$col <- colnames(C)[df.index$col]
    df.index$nrow <- unlist(tapply(df.index$row, df.index$row, function(x){1:length(x)}))

    df.index$coef <- C[which(C!=0)]
    df.index$coefname <- paste0(as.character(df.index$coef)," ")
    df.index[df.index$coef>0,"coefname"] <- paste0("+",df.index$coefname[df.index$coef>0])
    df.index[df.index$coefname == "+1 " & df.index$nrow == 1, "coefname"] <- ""
    df.index[df.index$coefname == "+1 ", "coefname"] <- " + "
            
    df.index[df.index$coefname == "-1 ", "coefname"] <- " - "

    df.index$rowname <- paste0(df.index$coefname,df.index$col)
    out <- unlist(tapply(df.index$rowname,df.index$row,paste,collapse=""))

    return(out)
}



##----------------------------------------------------------------------
### lTest.R ends here

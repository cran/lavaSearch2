### sCorrect-glht2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: apr 11 2023 (10:50) 
##           By: Brice Ozenne
##     Update #: 810
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## User 
### Code:

## * Documentation - glht2
#' @title General Linear Hypothesis Testing With Small Sample Correction
#' @description Test linear hypotheses on coefficients from a latent variable models with small sample corrections.
#' @name glht2
#' 
#' @param object,model a \code{lvmfit}, \code{lvmfit2}, or \code{mmm} object.
#' @param linfct [matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param robust [logical] should robust standard error be used? 
#' Otherwise rescale the influence function with the standard error obtained from the information matrix.
#' @param cluster  [integer vector] the grouping variable relative to which the observations are iid.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/code{FALSE}/code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... [logical] arguments passed to lower level methods.
#'
#' @details
#' Whenever the argument linfct is not a matrix, it is passed to the function \code{createContrast} to generate the contrast matrix and, if not specified, rhs. \cr \cr
#'
#' Since only one degree of freedom can be specify in a glht object and it must be an integer, the degree of freedom of the denominator of an F test simultaneously testing all hypotheses is retained, after rounding. \cr \cr
#'
#' Argument rhs and null are equivalent.
#' This redondance enable compatibility between \code{lava::compare}, \code{compare2}, \code{multcomp::glht}, and \code{glht2}.
#' @return A \code{glht} object.
#' 
#' @seealso
#' \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{estimate2}} to pre-compute quantities for the small sample correction.
#' 
#' @concept multiple comparisons
#'
#' @examples
#' library(multcomp)
#' 
#' ## Simulate data
#' mSim <- lvm(c(Y1,Y2,Y3)~ beta * eta, Z1 ~ E, Z2 ~ E, Age[40:5]~1)
#' latent(mSim) <- "eta"
#' set.seed(10)
#' n <- 1e2
#'
#' df.data <- lava::sim(mSim, n, latent = FALSE, p = c(beta = 1))
#'
#' #### Inference on a single model ####
#' e.lvm <- estimate(lvm(Y1~E), data = df.data)
#' summary(glht2(e.lvm, linfct = c("Y1~E + Y1","Y1")))
#' 
#' #### Inference on separate models ####
#' ## fit separate models
#' lvmX <- estimate(lvm(Z1 ~ E), data = df.data)
#' lvmY <- estimate(lvm(Z2 ~ E + Age), data = df.data)
#' lvmZ <- estimate(lvm(c(Y1,Y2,Y3) ~ eta, eta ~ E), 
#'                  data = df.data)
#'
#' #### create mmm object #### 
#' e.mmm <- mmm(X = lvmX, Y = lvmY, Z = lvmZ)
#'
#' #### create contrast matrix ####
#' resC <- createContrast(e.mmm, linfct = "E")
#'
#' #### adjust for multiple comparisons ####
#' e.glht2 <- glht2(e.mmm, linfct = c(X="E"), df = FALSE)
#' summary(e.glht2)
#'
#' @concept multiple comparison
#' @export
`glht2` <-
    function(object, ...) UseMethod("glht2")


## * glht2.lvmfit
#' @rdname glht2
#' @export
glht2.lvmfit <- function(object, linfct, rhs = NULL, robust = FALSE, cluster = NULL, ssc = lava.options()$ssc, df = lava.options()$df, ...){
    return(glht2(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...), linfct = linfct, rhs = rhs, robust = robust, cluster = cluster))

}

## * glht2.lvmfit2
#' @rdname glht2
#' @export
glht2.lvmfit2 <- function(object, linfct, rhs = NULL,
                          robust = FALSE, cluster = NULL,
                          ...){

    out <- compare2(object, linfct = linfct, rhs = rhs,
                    robust = robust, cluster = cluster,
                    as.lava = FALSE, F.test = FALSE, ...)
        
    return(out)
}


## * glht2.mmm
#' @rdname glht2
#' @export
glht2.mmm <- function (object, linfct, rhs = 0,
                       robust = FALSE, cluster = NULL,
                       ...){

    ## ** check the class of each model
    n.object <- length(object)
    name.object <- names(object)    
    if(is.null(name.object)){
        stop("Argument \'object\' must be named list. \n")
    }

    test.lvmfit <- sapply(object, inherits, what = "lvmfit")
    if(any(test.lvmfit == 0)){
        index.wrong <- which(test.lvmfit == 0)
        stop("Argument \'object\' must be a list of objects that inherits from lvmfit. \n",
             "Incorrect element(s): ",paste(index.wrong, collapse = " "),".\n")
    }
    test.lvmfit2 <- sapply(object, inherits, what = "lvmfit2")
    if(any(test.lvmfit2 == 0)){
        for(iO in which(test.lvmfit2==0)){
            object[[iO]] <- estimate2(object[[iO]], dVcov.robust = robust, ...)
        }
    }
    
    ## ** define the contrast matrix
    out <- list()
    if (is.character(linfct)){
        resC <- createContrast(object, linfct = linfct, rowname.rhs = FALSE)
        linfct <- resC$contrast
        ls.contrast <- resC$mlf
        if("rhs" %in% names(match.call()) == FALSE){
            rhs <- resC$null
        }
    }else if(is.matrix(linfct)){
        ls.contrast <- lapply(name.object, function(x){ ## x <- name.object[2]
            iColnames <- grep(paste0("^",x,": "), colnames(linfct), value = FALSE, fixed = FALSE)
            iRownames <- rowSums(linfct[,iColnames]!=0)>0
            linfct[iRownames, iColnames,drop=FALSE]            
        })
        names(ls.contrast) <- name.object
        contrast <- linfct
        if("rhs" %in% names(match.call()) == FALSE){ ## left rhs to default value
            rhs <- rep(0, NROW(contrast))
        }else if(length(rhs)!=NROW(contrast)){
            stop("mismatch between the dimensions of argument \'rhs\' and argument \'contrast\' \n")
        }
    }else{
        stop("Argument \'linfct\' must be a matrix or a vector of characters. \n",
             "Consider using  out <- createContrast(...) and pass out$contrast to linfct. \n")
    }

    ## ** check whether it is possible to compute df
    ## i.e. are linear hypothesis model specific?
    test.df <- all(unlist(lapply(object, function(iModel){iModel$sCorrect$df == "satterthwaite"})))
    if(test.df){
        n.hypo <- NROW(linfct)
        ls.modelPerTest <- lapply(1:n.hypo, function(iHypo){ ## iHypo <- 1
            iContrast <- linfct[iHypo,]
            iNames <- names(iContrast)[abs(iContrast)>0]
            iModels <- unlist(lapply(strsplit(iNames, split = ":"),"[[",1))
            return(length(unique(iModels)))
        })
        
        if(any(unlist(ls.modelPerTest)>1)){
            stop("Cannot compute the degrees of freedom for tests performed across several models \n",
                 "Consider setting the argument \'df\' to FALSE \n")
        }    
    }

    ## ** Total number of observations
    if(!is.null(cluster)){
        ls.cluster <- lapply(object, function(iO){extractData(iO, rm.na = FALSE)[[cluster]]})
        Ucluster <- unique(unlist(ls.cluster))
        n.cluster <- length(Ucluster)
    }
    
    ## ** Extract influence functions from all models
    ls.res <- lapply(1:n.object, function(iM){ ## iM <- 1

        ## *** Pre-compute quantities
        if(!inherits(object[[iM]],"lvmfit2")){
            object[[iM]] <- estimate2(object[[iM]], ...)
        }
        out$param <- coef(object[[iM]], as.lava = FALSE)
        name.param <- names(out$param)
        name.object.param <- paste0(name.object[iM],": ",name.param)
        out$param <- stats::setNames(out$param, name.object.param)
        
        ## *** Compute df for each test
        if(!is.na(object[[iM]]$sCorrect$df)){
            ## here null does not matter since we only extract the degrees of freedom
            iContrast <- ls.contrast[[iM]]
            colnames(iContrast) <- name.param
            iWald <- compare2(object[[iM]], linfct = iContrast, as.lava = FALSE, F.test = FALSE)
            out$df <- iWald$df
        }else{
            out$df <- Inf
        }
        ## *** get iid decomposition
        iid.tempo <- iid(object[[iM]], robust = robust, cluster = cluster, as.lava = FALSE)
        if(!is.null(cluster)){
            out$iid <- matrix(NA, nrow = n.cluster, ncol = length(name.param),
                              dimnames = list(Ucluster, name.param))
            out$iid[attr(iid.tempo,"cluster"),] <- iid.tempo
        }else{
            out$iid <- iid.tempo
        }
        colnames(out$iid) <- name.object.param

        ## *** get se
        if(robust){
            out$se <- sqrt(diag(crossprod(iid.tempo)))
        }else{
            out$se <- sqrt(diag(vcov(object[[iM]], as.lava = FALSE)))
        }
        return(out)
        
    })
    seq.df <- unlist(lapply(ls.res,"[[","df"))
    seq.param <- unlist(lapply(ls.res,"[[","param"))

    if(test.df){
        df.global <- round(stats::median(seq.df), digits = 0)
    }else{
        df.global <- 0
    }
    ls.iid <- lapply(ls.res,"[[","iid")
    ls.se <- lapply(ls.res,"[[","se")
    n.obs <- unique(unlist(lapply(ls.iid, NROW)))
    if(length(n.obs)>1){
        stop("Mismatch between the number of observations in the iid \n",
             "Likely to be due to the presence of missing values \n",
             "Consider specifying the \'cluster\' argument \n")
    }
    M.iid <- do.call(cbind,ls.iid)
    diag.se <- diag(do.call(c,ls.se))
    if(any(is.na(M.iid))){
       M.iid[is.na(M.iid)] <- 0
    }
    vcov.object <- diag.se %*% stats::cov2cor(crossprod(M.iid)) %*% diag.se ## same as multcomp:::vcov.mmm
    dimnames(vcov.object) <- list(colnames(M.iid), colnames(M.iid))
    
    ## ** sanity check
    name.param <- names(seq.param)
    if(!identical(colnames(linfct),name.param)){
        stop("Column names of the contrast matrix does not match the one of the coefficients \n")
    }
    if(!identical(colnames(vcov.object),name.param)){
        stop("Column names of the variance covariance matrix does not match the one of the coefficients \n")
    }
    if(!identical(rownames(vcov.object),name.param)){
        stop("Rownames names of the variance covariance matrix does not match the one of the coefficients \n")
    
    }

    ## ** convert to the appropriate format    
    out <- list(model = object,
                linfct = linfct,
                rhs = unname(rhs),
                coef = seq.param,
                vcov = vcov.object,
                df = df.global,
                alternative = "two.sided",
                type = NULL,
                robust = robust)
    class(out) <- c("glht2","glht")
        
    ### ** export
    return(out)    
}


## * glht.lvmfit2
#' @rdname glht2
#' @export
glht.lvmfit2 <- function(model, linfct, rhs = NULL,
                         robust = FALSE, cluster = NULL,
                         ...){

    out <- compare2(model, linfct = linfct, rhs = rhs,
                    robust = robust, cluster = cluster,
                    as.lava = FALSE, F.test = FALSE, ...)
        
    return(out)
}

## * .calcClosure
.calcClosure <- function(name, estimate, covariance, type, df){

    n.hypo <- length(name)
    correlation <- stats::cov2cor(covariance)

    ## ** create all possible hypotheses
    ls.closure <- lapply(n.hypo:1, function(iNtest){ ## iNtest <- 1  
        iList <- list(M = utils::combn(name, m = iNtest))
        iList$vec <- apply(iList$M, 2, paste, collapse = ",")
        return(iList)
    })

    ## ** compute all p.values
    for(iLevel in 1:length(ls.closure)){ ## iLevel <- 1
        ls.closure[[iLevel]]$test <- t(apply(ls.closure[[iLevel]]$M, 2, function(iHypo){
            index <- which(name %in% iHypo)
            if(type == "chisq"){
                return(.ChisqTest(estimate[index], covariance = covariance[index,index,drop=FALSE], df = df))
            }else if(type == "max"){
                return(.tTest(estimate[index],
                              covariance = covariance[index,index,drop=FALSE],
                              correlation = correlation[index,index,drop=FALSE], df = df))
            }
        }))
        rownames(ls.closure[[iLevel]]$test) <- ls.closure[[iLevel]]$vec
    }
    
    ## ** find all hypotheses in the closure related to an individual hypothesis
    ls.hypo <- vector(mode = "list", length = n.hypo)
    for(iHypo in 1:n.hypo){ ## iHypo <- 1
        ls.hypo[[iHypo]] <- do.call(rbind,lapply(ls.closure, function(iClosure){ ## iClosure <- 1
            iIndex <- which(colSums(iClosure$M==name[iHypo])>0)
            data.frame(hypothesis = iClosure$vec[iIndex],
                       statistic = as.double(iClosure$test[iIndex,"statistic"]),
                       p.value = as.double(iClosure$test[iIndex,"p.value"]))
        }))
    }
    names(ls.hypo) <- name
        
    ## ** adjusted p.values
    vec.p.value <- unlist(lapply(ls.hypo, function(x){max(x$p.value)}))
    return(list(closure = ls.closure,
                test = ls.hypo,
                p.value = vec.p.value))
    
}

## * .tTest
.tTest <- function(estimate, covariance, correlation, df, ...){
    df1 <- length(estimate)
    statistic <- max(abs(estimate/sqrt(diag(covariance))))
    if(is.null(df)){
        distribution <-  "gaussian"
    }else{
        distribution <- "student"
    }
    p.value <- .calcPmaxIntegration(statistic, p = df1, Sigma = correlation, df = df,
                                    distribution = distribution)
    return(c("statistic" = statistic,
             "p.value" = p.value))
}

## * .ChisqTest
.ChisqTest <- function(estimate, covariance, df, ...){
    df1 <- length(estimate)
    ## q * statistic ~ chisq or fisher
    statistic <- as.double(matrix(estimate, nrow = 1) %*% solve(covariance) %*% matrix(estimate, ncol = 1)) / df1
    if(!is.null(df)){
        return(c("statistic" = statistic,
                 "p.value" = 1-stats::pf(statistic, df1 = df1, df2 = df)))
    }else{
        return(c("statistic" = statistic,
                 "p.value" = 1-stats::pchisq(statistic, df = df1)))
        
    }
}


 

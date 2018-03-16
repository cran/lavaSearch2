### mlf2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: mar 12 2018 (16:59) 
##           By: Brice Ozenne
##     Update #: 358
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## User 
### Code:


## * estfun.lvmfit
#' @title Extract Empirical Estimating Functions (lvmfit Object)
#' @description Extract the empirical estimating functions of a lvmfit object.
#' This function is for internal use but need to be public to enable its use by \code{multcomp::glht}.
#' 
#' @param x an \code{lvmfit} object.
#' @param ... arguments passed to methods.
#'
#' @details This function enables to use the \code{glht} function with lvmfit object.
#' Otherwise when calling \code{multcomp:::vcov.mmm} then \code{sandwich::sandwich} and then \code{sandwich::meat}, \code{sandwich::meat} will complain that \code{estfun} is not defined for \code{lvmfit} objects.
#'
#' @examples
#' library(multcomp)
#' 
#' #### generative model ####
#' mSim <- lvm(X ~ Age + 0.5*Treatment,
#'             Y ~ Gender + 0.25*Treatment,
#'             c(Z1,Z2,Z3) ~ eta, eta ~ 0.75*treatment,
#'             Age[40:5]~1)
#' latent(mSim) <- ~eta
#' categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
#' categorical(mSim, labels = c("male","female")) <- ~Gender
#'
#' #### simulate data ####
#' n <- 5e1
#' set.seed(10)
#' df.data <- sim(mSim, n = n, latent = FALSE)
#'
#' #### fit separate models ####
#' lmX <- lm(X ~ Age + Treatment, data = df.data)
#' lvmY <- estimate(lvm(Y ~ Gender + Treatment), data = df.data)
#' lvmZ <- estimate(lvm(c(Z1,Z2,Z3) ~ eta, eta ~ Treatment), 
#'                  data = df.data)
#'
#' #### create mmm object #### 
#' e.mmm <- mmm(X = lmX, Y = lvmY, Z = lvmZ)
#'
#' #### create contrast matrix ####
#' resC <- createContrast(e.mmm, var.test = "Treatment", add.variance = FALSE)
#'
#' #### adjust for multiple comparisons ####
#' e.glht <- glht(e.mmm, linfct = resC$mlf)
#' summary(e.glht)
#' 
#' @method estfun lvmfit
#' @concept multiple comparison
#' @export
estfun.lvmfit <- function(x, ...){
    U <- lava::score(x, indiv = TRUE)
    return(U)
}


## * Documentation - glht2
#' @title General Linear Hypothesis
#' @description Test general linear hypotheses and across latent variable models with small sample corrections.
#' @name glht2
#' 
#' @param model a \code{lvmfit} or \code{mmm} object.
#' The \code{mmm} object can only contain lm/gls/lme/lvmfit objects.
#' @param linfct [matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param robust [logical] should robust standard error be used? 
#' Otherwise rescale the influence function with the standard error obtained from the information matrix.
#' @param ... [internal] Only used by the generic method.
#'
#' @details
#' Whenever the argument linfct is not a matrix, it is passed to the function \code{createContrast} to generate the contrast matrix and, if not specified, rhs. \cr \cr
#'
#' Since only one degree of freedom can be specify in a glht object and it must be an integer, the degree of freedom of the denominator of an F test simultaneously testing all hypotheses is retained, after rounding.
#'
#' @return A \code{glht} object.
#' 
#' @seealso
#' \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{sCorrect}} to pre-compute quantities for the small sample correction.
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
#' df.data <- sim(mSim, n, latent = FALSE, p = c(beta = 1))
#'
#' #### Inference on a single model ####
#' e.lvm <- estimate(lvm(Y1~E), data = df.data)
#' summary(glht2(e.lvm, linfct = c("Y1~E + Y1","Y1")))
#' 
#' #### Inference on separate models ####
#' ## fit separate models
#' lmX <- lm(Z1 ~ E, data = df.data)
#' lvmY <- estimate(lvm(Z2 ~ E + Age), data = df.data)
#' lvmZ <- estimate(lvm(c(Y1,Y2,Y3) ~ eta, eta ~ E), 
#'                  data = df.data)
#'
#' #### create mmm object #### 
#' e.mmm <- mmm(X = lmX, Y = lvmY, Z = lvmZ)
#'
#' #### create contrast matrix ####
#' resC <- createContrast(e.mmm, var.test = "E", add.variance = TRUE)
#'
#' #### adjust for multiple comparisons ####
#' e.glht2 <- glht2(e.mmm, linfct = resC$contrast)
#' summary(e.glht2)
#'
#' @concept multiple comparison
#' @export
`glht2` <-
  function(model, ...) UseMethod("glht2")

## * glht2.lvmfit
#' @rdname glht2
#' @export
glht2.lvmfit <- function(model, linfct, rhs = 0,
                         bias.correct = TRUE, robust = FALSE, ...){

    ### ** define contrast matrix
    if(!is.matrix(linfct)){
        resC <- createContrast(model, par = linfct)
        linfct <- resC$contrast
        if("rhs" %in% names(match.call()) == FALSE){
            rhs <- resC$null
        }
    }

### ** pre-compute quantities for the small sample correction
    if(!inherits(model,"lvmfit2")){
        sCorrect(model, score = robust) <- bias.correct
    }

### ** Wald test with small sample correction
    name.param <- colnames(linfct)
    n.param <- NCOL(linfct)
    n.hypo <- NROW(linfct)

    resWald <- compare2(model, contrast = linfct, null = rhs, as.lava = FALSE)
    ## update name according to multcomp, i.e. without second member
    rownames(linfct) <- .contrast2name(linfct, null = NULL) 

    ### ** Global degree of freedom
    df.global <- round(resWald["global","df"], digits = 0)
    
    ### ** compute variance-covariance matrix
    if(robust){
        vcov.model <- crossprod(iid2(model))
    }else{
        vcov.model <- model$sCorrect$vcov.param
    }

### ** convert to the appropriate format
    out <- list(model = model,
                linfct = linfct,
                rhs = unname(rhs),
                coef = coef(model),
                vcov = vcov.model,
                df = df.global,
                alternative = "two.sided",
                type = NULL)
    class(out) <- "glht"
        
    ### ** export
    return(out)
}


## * glht2.mmm
#' @rdname glht2
#' @export
glht2.mmm <- function (model, linfct, rhs = 0, bias.correct = TRUE, robust = FALSE, ...){
    ### ** check the class of each model
    n.model <- length(model)
    name.model <- names(model)    
    if(is.null(name.model)){
        stop("Argument \'model\' must be named list. \n")
    }
    
    test.lm <- sapply(model, inherits, what = "lm")
    test.gls <- sapply(model, inherits, what = "gls")
    test.lme <- sapply(model, inherits, what = "lme")
    test.lvmfit <- sapply(model, inherits, what = "lvmfit")
    if(any(test.lm + test.gls + test.lme + test.lvmfit == 0)){
        index.wrong <- which(test.lm + test.gls + test.lme + test.lvmfit == 0)
        stop("Argument \'model\' must be a list of objects that inherits from lm/gls/lme/lvmfit. \n",
             "Incorrect element(s): ",paste(index.wrong, collapse = " "),".\n")
    }

    ### ** define the contrast matrix
    out <- list()
    if (is.character(linfct)){
        resC <- createContrast(model, par = linfct, add.variance = TRUE)
        contrast <- resC$contrast
        ls.contrast <- resC$ls.contrast
        if("rhs" %in% names(match.call()) == FALSE){
            rhs <- resC$null
        }
    }else if(is.matrix(linfct)){
        
        ls.contrast <- lapply(name.model, function(x){ ## x <- name.model[2]
            
            iRownames <- grep(paste0(x,": "), rownames(linfct), value = FALSE, fixed = TRUE)
            iColnames <- grep(paste0(x,": "), colnames(linfct), value = FALSE, fixed = TRUE)
            linfct[iRownames, iColnames,drop=FALSE]            
        })
        names(ls.contrast) <- name.model
        contrast <- linfct
    }else{
        stop("Argument \'linfct\' must be a matrix or a vector of characters. \n",
             "Consider using  out <- createContrast(...) and pass out$contrast to linfct. \n")
    }

    ## ** Extract influence functions from all models    
    ls.res <- lapply(1:n.model, function(iM){ ## iM <- 1

### *** Pre-compute quantities
        if(!inherits(model[[iM]],"lm2") && !inherits(model[[iM]],"gls2") && !inherits(model[[iM]],"lme2") && !inherits(model[[iM]],"lvmfit2")){
            sCorrect(model[[iM]], score = TRUE) <- bias.correct
        }
        out$param <- model[[iM]]$sCorrect$param
        name.param <- names(out$param)
        name.object.param <- paste0(name.model[iM],": ",name.param)
        out$param <- setNames(out$param, name.object.param)
        
### *** Compute df for each test
        ## here null does not matter since we only extract the degrees of freedom
        iContrast <- ls.contrast[[iM]]
        colnames(iContrast) <- name.param
        
        iWald <- compare2(model[[iM]], contrast = iContrast, as.lava = FALSE)
        out$df <- iWald[1:(NROW(iWald)-1),"df"]

### *** get iid decomposition
        out$iid <- iid2(model[[iM]], robust = robust)
        colnames(out$iid) <- name.object.param
            
        return(out)
        
    })
    seq.df <- unlist(lapply(ls.res,"[[","df"))
    seq.param <- unlist(lapply(ls.res,"[[","param"))
    df.global <- round(stats::median(seq.df), digits = 0)
    vcov.model <- crossprod(do.call(cbind,lapply(ls.res,"[[","iid")))

    ### ** convert to the appropriate format
    out <- list(model = model,
                linfct = linfct,
                rhs = unname(rhs),
                coef = seq.param,
                vcov = vcov.model,
                df = df.global,
                alternative = "two.sided",
                type = NULL)
    class(out) <- "glht"
        
    ### ** export
    return(out)    
}


### createContrast.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 31 2018 (12:05) 
## Version: 
## Last-Updated: mar 12 2018 (18:01) 
##           By: Brice Ozenne
##     Update #: 200
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - createContrast
#' @title Create Contrast matrix
#' @description Returns a contrast matrix corresponding an object.
#' The contrast matrix will contains the hypotheses in rows and the model coefficients in columns.
#' @name createContrast
#' 
#' @param object a \code{ls.lvmfit} object.
#' @param par [vector of characters] expression defining the linear hypotheses to be tested. See the examples section. 
#' @param add.variance [logical] should the variance coefficients be considered as model coefficients?
#' Required for lm, gls, and lme models.
#' @param var.test [character] a regular expression that is used to identify the coefficients to be tested using \code{grep}. Each coefficient will be tested in a separate hypothesis. When this argument is used, the argument \code{par} is disregarded.
#' @param name.param [internal] the names of all the model coefficients.
#' @param add.rowname [internal] should a name be defined for each hypothesis.
#' @param rowname.rhs should the right hand side of the null hypothesis be added to the name.
#' @param ... [internal] Only used by the generic method.
#'
#' @details
#' One can initialize an empty contrast matrix setting the argument\code{par} to \code{character(0)}. \cr \cr
#'
#' When using \code{multcomp::glht} one should set the argument \code{add.variance} to \code{FALSE}. \cr
#' When using \code{lavaSearch2::glht2} one should set the argument \code{add.variance} to \code{TRUE}.
#' 
#' @return A list containing
#' \itemize{
#' \item{contrast} [matrix] a contrast matrix corresponding to the left hand side of the linear hypotheses.
#' \item{null} [vector] the right hand side of the linear hypotheses.
#' \item{Q} [integer] the rank of the contrast matrix.
#' \item{ls.contrast} [list, optional] the contrast matrix corresponding to each submodel.
#' Only present when the \code{argument} object is a list of models.
#' }
#' @examples
#' ## Simulate data
#' mSim <- lvm(X ~ Age + Treatment,
#'             Y ~ Gender + Treatment,
#'             c(Z1,Z2,Z3) ~ eta, eta ~ treatment,
#'             Age[40:5]~1)
#' latent(mSim) <- ~eta
#' categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
#' categorical(mSim, labels = c("male","female")) <- ~Gender
#' n <- 1e2
#' set.seed(10)
#' df.data <- sim(mSim,n)
#'
#' ## Estimate separate models
#' lmX <- estimate(lvm(X ~ -1 + Age + Treatment), data = df.data)
#' lmY <- estimate(lvm(Y ~ -1 + Gender + Treatment), data = df.data)
#' lvmZ <- estimate(lvm(c(Z1,Z2,Z3) ~ -1 + 1*eta, eta ~ -1 + Treatment), 
#'                  data = df.data)
#'
#' ## Contrast matrix for a given model
#' createContrast(lmX, par = "X~Age")
#' createContrast(lmX, par = c("X~Age=0","X~Age+5*X~TreatmentSSRI=0"))
#' createContrast(lmX, par = character(0))
#'
#' ## Contrast matrix for the join model
#' ls.lvm <- list(X = lmX, Y = lmY, Z = lvmZ)
#' createContrast(ls.lvm, var.test = "Treatment")
#' createContrast(ls.lvm, par = character(0))
#'
#' @concept small sample inference
#' 
#' @export
`createContrast` <-
    function(object, ...) UseMethod("createContrast")

## * createContrast.character
#' @rdname createContrast
#' @export
createContrast.character <- function(object, name.param,
                                     add.rowname = TRUE, rowname.rhs = TRUE,
                                     ...){

    n.param <- length(name.param)
    
    n.hypo <- length(object)
    if(any(nchar(object)==0)){
        stop("Argument contains empty character string(s) instead of an expression involving the model mean coefficients \n")
    }
    null <- rep(NA, n.hypo)
    contrast <- matrix(0, nrow = n.hypo, ncol = n.param,
                       dimnames = list(NULL,name.param))

    if(n.hypo>0){
        for(iH in 1:n.hypo){ # iH <- 1
            iTempo.eq <- strsplit(object[iH], split = "=", fixed = TRUE)[[1]]
            if(length(iTempo.eq)==1){ ## set null to 0 when second side of the equation is missing
                iTempo.eq <- c(iTempo.eq,"0")
            }
            null[iH] <- as.numeric(trim(iTempo.eq[2]))
            iRh.plus <- strsplit(iTempo.eq[[1]], split = "+", fixed = TRUE)[[1]]
            iRh <- trim(unlist(sapply(iRh.plus, strsplit, split = "-", fixed = TRUE)))
            iRh <- iRh[iRh!="",drop=FALSE]
                            
            ls.iRh <- lapply(strsplit(iRh, split = "*", fixed = TRUE), trim)
                    
            iN.tempo <- length(ls.iRh)
                    
            for(iCoef in 1:iN.tempo){ # iCoef <- 1

                if(length(ls.iRh[[iCoef]])==1){
                    iFactor <- 1
                    iName <- ls.iRh[[iCoef]][1]                
                }else{
                    iFactor <- as.numeric(ls.iRh[[iCoef]][1])
                    iName <- ls.iRh[[iCoef]][2]
                }
            
                if(iName %in% name.param == FALSE){
                    txt.message <- paste0("unknown coefficient ",iName," in hypothesis ",iH,"\n")
                    possibleMatch <- pmatch(iName, table = name.param)
                    if(identical(possibleMatch, as.integer(NA))){
                        possibleMatch <- grep(iName, name.param, fixed = TRUE, value = TRUE)
                    }
                    if(!identical(possibleMatch, as.integer(NA))){
                        txt.message <- c(txt.message,
                                         paste0("candidates: \"",paste(possibleMatch, collapse = "\" \""),"\"\n"))
                    }
                    stop(txt.message)                    
                }

                test.sign <- length(grep("-",strsplit(iRh[iCoef], split = iName)[[1]][1]))>0
                contrast[iH,iName] <- c(1,-1)[test.sign+1] * iFactor
            }
        }
    
        if(add.rowname){
            name.hypo <- .contrast2name(contrast, null = if(rowname.rhs){null}else{NULL})
            rownames(contrast) <- name.hypo
            null <- setNames(null, name.hypo)
        }
    }
    
    return(list(contrast = contrast,
                null = null,
                Q = n.hypo))
}

## * createContrast.lm
#' @rdname createContrast
#' @export
createContrast.lm <- function(object, par, add.variance, ...){

    if(!identical(class(par),"character")){
        stop("Argument \'par\' must be a character \n")
    }    
    name.coef <- names(coef(object))
    if(add.variance){
        if(any("sigma2" %in% name.coef)){
            stop("createContrast does not work when one of the coefficients is named \"sigma2\" \n")
        }
        name.coef <- c(name.coef,"sigma2")
    }
    
    out <- createContrast(par, name.param = name.coef, ...)
    return(out)
    
}

## * createContrast.gls
#' @rdname createContrast
#' @export
createContrast.gls <- function(object, par, add.variance, ...){

    if(!identical(class(par),"character")){
        stop("Argument \'par\' must be a character \n")
    }
    if(add.variance){
        name.coef <- names(.coef2(object))
    }else{
        name.coef <- names(coef(object))
    }
    out <- createContrast(par, name.param = name.coef, ...)
    return(out)
    
}

## * createContrast.lme
#' @rdname createContrast
#' @export
createContrast.lme <- createContrast.lm

## * createContrast.lvmfit
#' @rdname createContrast
#' @export
createContrast.lvmfit <- function(object, par, ...){

    if(!identical(class(par),"character")){
        stop("Argument \'par\' must be a character \n")
    }
    out <- createContrast(par, name.param = names(coef(object)), ...)
    return(out)
    
}

## * createContrast.list
#' @rdname createContrast
#' @export
createContrast.list <- function(object, par = NULL, add.variance, var.test = NULL, 
                                ...){

    ## ** find the names of the coefficients
    name.model <- names(object)
    
    ls.coefname <- lapply(name.model, function(iModel){ ## list by model
        iResC <- createContrast(object[[iModel]], par = character(0), add.variance = add.variance)
        return(colnames(iResC$contrast))
    })
    names(ls.coefname) <- name.model

    ls.object.coefname <- lapply(name.model, function(iModel){ ## list by model with model name
        paste0(iModel,": ", ls.coefname[[iModel]])
    })    
    names(ls.object.coefname) <- name.model
    
    object.coefname <- unname(unlist(ls.object.coefname)) ## vector
    n.coef <- length(object.coefname)
    
    ## ** normalize arguments
    if(!is.null(var.test)){
        if(!is.null(par)){
            stop("Argument \'var.test\' cannot be specified when argument \'par\' is specified \n")
        }else{
            if(length(var.test)!=1){
                stop("Argument \'var.test\' must have length 1 \n")
            }
            par <- grep(var.test, object.coefname, value = TRUE)
        }
    }

    ## ** create full contrast matrix
    out <- createContrast(par, name.param = object.coefname)
    if(any(out$null!=0)){
        warning("glht ignores the \'rhs\' argument when dealing with a multiple models \n")
    }

    ## ** create contrast matrix relative to each model
    out$mlf <- lapply(name.model, function(iModel){ ## iModel <- name.model[1]        
        ## only keep columns corresponding to coefficients belonging the the current model
        iContrast <- out$contrast[,ls.object.coefname[[iModel]],drop=FALSE]

        ## update name by removing the name of the model
        colnames(iContrast) <- ls.coefname[[iModel]]

        ## remove lines in the contrast matrix containing only 0
        if(NROW(iContrast)>0){
          index.n0 <- which(rowSums(iContrast!=0)!=0)
          return(iContrast[index.n0,,drop=FALSE])
        }else{
          return(iContrast)
        }
    })
    names(out$mlf) <- name.model    
    class(out$mlf) <- "mlf"
    
    ## remove right hand side from the names (like in multicomp)
    if(length(par)>0){
        rownames(out$contrast) <- .contrast2name(out$contrast, null = NULL)
        out$mlf <- lapply(out$mlf, function(x){ ## x <- name.model[1]
            rownames(x) <- .contrast2name(x, null = NULL)
            return(x)
        })
        class(out$mlf) <- "mlf"
        names(out$null) <- rownames(out$contrast)
    }
   
    ## ** export
    return(out)    
}

## * createContrast.mmm
#' @rdname createContrast
#' @export
createContrast.mmm <- createContrast.list

## * .contrast2name
#' @title Create Rownames for a Contrast Matrix
#' @description Create rownames for a contrast matrix using the coefficients and the names of the coefficients. The rownames will be [value * name] == null, e.g. [beta + 4*alpha] = 0.
#' @name contrast2name
#'
#' @param contrast [matrix] a contrast matrix defining the left hand side of the linear hypotheses to be tested.
#' @param null [vector, optional] the right hand side of the linear hypotheses to be tested.
#'
#' @details When argument \code{NULL} is null then the rownames will not be put into brackets and the right hand side will not be added to the name.
#'
#' @return a character vector.
#' 
#' @keywords internal
.contrast2name <- function(contrast, null = NULL){

    contrast.names <- colnames(contrast)
    
    df.index <- as.data.frame(which(contrast != 0, arr.ind = TRUE))
    df.index$col <- contrast.names[df.index$col]
    index.order <- order(df.index$row,match(df.index$col,contrast.names))
    df.index <- df.index[index.order,]
    df.index$nrow <- unlist(tapply(df.index$row, df.index$row, function(x){1:length(x)}))

    ## find coef value to each  coefficient
    df.index$coef <- contrast[which(contrast!=0)][index.order]
    df.index$coefname <- as.character(df.index$coef)

    ## add positive sign
    df.index[df.index$coef>0,"coefname"] <- paste0("+",df.index$coefname[df.index$coef>0])
    
    ## add multiplicative symbol
    index.Npm1 <- which(df.index$coefname %in% c("","+1","-1") == FALSE)
    df.index[index.Npm1, "coefname"] <- paste0(df.index[index.Npm1, "coefname"],"*")

    ## simplify
    df.index[df.index$coefname == "+1" & df.index$nrow == 1, "coefname"] <- ""
    df.index[df.index$coefname == "+1", "coefname"] <- "+"            
    df.index[df.index$coefname == "-1", "coefname"] <- "-"
    df.index[df.index$nrow == 1, "coefname"] <- gsub("+","",df.index[df.index$nrow == 1, "coefname"], fixed = TRUE)

    ## add space between coefficients
    df.index$coefname <- gsub("+"," + ",df.index$coefname, fixed = TRUE)
    df.index$coefname <- gsub("-"," - ",df.index$coefname, fixed = TRUE)
    df.index[df.index$coefname == " - " & df.index$nrow == 1, "coefname"] <- "- "

    ## paste together value and coefficient names
    df.index$rowname <- paste0(df.index$coefname,df.index$col)
    ## paste together names from the same hypothesis
    out <- unlist(tapply(df.index$rowname,df.index$row,paste,collapse=""))

    ## add right hand side
    if(!is.null(null)){
        out <- paste0("[",out,"] = ",null)
        
    }

    return(as.character(out))
}


##----------------------------------------------------------------------
### createContrast.R ends here

### var2dummy.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 22 2017 (16:03) 
## Version: 
## last-updated: jan 16 2018 (09:06) 
##           By: Brice Ozenne
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Convert Variable Names to Dummy Variables Names.
#' @description When dealing with categorical variables, the \code{estimate} function convert the categorical variables into dummy variables.
#' This function convert a set of variable names to their corresponding name in the model with dummy variables
#' @name var2dummy
#' 
#' @param x a latent variable model.
#' @param var the variable to be transformed.
#' @param data dataset according to which the model should be updated.
#' @param rm.first.factor should the first level of each categorical variable be ignored?
#' @param ... additional arguments to be passed to lower levels functions.
#' 
#' @examples
#' library(lava)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u ~ X1+X2
#' var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("M","F","MF")) <- ~X1
#' var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("1","2","3")) <- ~X2
#' var2dummy(m, var = c("X1","X2"))
#' @export
`var2dummy` <-
  function(x,...) UseMethod("var2dummy")

#' @rdname var2dummy
#' @export
var2dummy.list <- function(x, var, rm.first.factor = TRUE, ...){

    var <- stats::setNames(var,var)
    ## convertion to dummy variable name for categorical variables
    factor.var <- names(x$x$attributes$labels)
    
    if(!is.null(var) && any(var %in% factor.var)){
        subvar <- var[var %in% factor.var]
        for(iFactor in subvar){ # iFactor <- "X1"
            newvar <- paste0(iFactor,x$x$attributes$labels[[iFactor]])
            if(rm.first.factor){newvar <- newvar[-1]}
            newvar <- stats::setNames(newvar, rep(iFactor, length(newvar)))
            var <- c(var[names(var)!=iFactor],newvar)            
        }
    }
    return(var)
}

#' @rdname var2dummy
#' @export
var2dummy.lvm <- function(x, data = NULL, ...){

    if(is.null(data)){
        data <- lava::sim(x,1)
    }

    x2 <- lava_categorical2dummy(x, data)

    ## recover attributes for models not defined using categorical
    obsvars <- setdiff(vars(x2$x),latent(x))
    if(any(obsvars %in% names(data) == FALSE)){
        missing.vars <- obsvars[obsvars %in% names(data) == FALSE]
        test.num <- sapply(1:NCOL(data), function(col){is.numeric(data[[col]])})
        possible.match <- names(data)[test.num==FALSE]
        n.possible.match <- length(possible.match)
        
        ls.labels <- list()
        for(iMatch in 1:n.possible.match){ # iMatch <- 1
            iVar <- possible.match[iMatch]
            iLabel <- levels(as.factor(data[[iVar]]))
            if(any(paste0(iVar,iLabel) %in% missing.vars)){
                ls.labels[[iVar]] <- iLabel
            }
        }

        x2$x$attributes$labels <- ls.labels
    }
    
    res <- var2dummy(x2, ...)
    return(res)
}

#----------------------------------------------------------------------
### var2dummy.R ends here

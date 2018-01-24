## * Documentation - initVarLink
#' @title Normalize var1 and var2
#' @name initVarLink
#' @description Convert var1 and var2 from formula or covariance to character.
#' 
#' @param var1 a character indicating the endogeneous variable or a formula
#' @param var2 an optional character indicating the exogeneous variable
#' @param repVar1 should var1 be duplicated to match var2 length. Only active if format = "list".
#' @param format should the name of the variable be return (format = "list"), a vector of character formula ("txt.formula") or a list of formula ("formula")
#' @param Slink the symbol for regression link
#' @param Scov the symbol for covariance link
#' @param ... argument to be passed to lower level functions
#' 
#' @details See test file test/testthat/test-Reduce.R for examples
#'
#' @examples
#' initVarLink(y ~ x1)
#' initVarLink("y ~ x1")
#' initVarLink(y ~ x1 + x2)
#' initVarLink("y ~ x1 + x2")
#' initVarLink(y ~ x1 + x2, repVar1 = TRUE)
#' initVarLink(y ~ x1 + x2, repVar1 = TRUE, format = "formula")
#' initVarLink(y ~ x1 + x2, repVar1 = TRUE, format = "txt.formula")
#' initVarLink("y", "x1", format = "formula")
#'
#' initVarLink("y ~ x1:0|1")
#'
#' initVarLinks(y ~ x1)
#' initVarLinks("y ~ x1")
#' initVarLinks(c("y ~ x1","y~ x2"))
#' initVarLinks(c(y ~ x1,y ~ x2))
#' initVarLinks(c("y ~ x1","y ~ x2"), format = "formula")
#' initVarLinks(c(y ~ x1,y ~ x2), format = "formula")
#' initVarLinks(c("y ~ x1","y~ x2"), format = "txt.formula")
#' initVarLinks(c(y ~ x1,y ~ x2), format = "txt.formula")

## * initVarLink
#' @rdname initVarLink
#' @export
initVarLink <- function(var1, var2, repVar1 = FALSE, format = "list",
                         Slink = c(lava.options()$symbols[1],"~"),
                         Scov = lava.options()$symbols[2]){

    format <- match.arg(format, c("list","txt.formula","formula"))
    test.formula <- (class(var1) == "formula")
    test.covariance <- sapply(Scov,grepl,x=var1,fixed=TRUE)
    test.regression <- sapply(Slink,grepl,x=var1,fixed=TRUE)

    if(missing(var2)){
        if(test.formula){
            var2 <- selectRegressor(var1, type = "vars")
            var1 <- selectResponse(var1, type = "vars")
            sep <- if(format == "formula"){"~"}else{Slink}
        }else if(any(test.covariance)){ ## covariance
            Scov <- Scov[test.covariance][1]
            varSplit <- strsplit(var1, split = Scov)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Scov}
        } else if(any(test.regression)){ ## regression
            Slink <- Slink[test.regression][1]
            varSplit <- strsplit(var1, split = Slink)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Slink}
        } else {
            var1 <- var1
            var2 <- NA
        }
    }else{

        if(!is.character(var1) || !is.character(var2)){
            stop("\'var1\' and \'var2\' must be characters when both are specified \n")
        }
        sep <- if(format == "formula"){"~"}else{Slink}        
    }
  
  
#### convert to format
    if(format == "formula"){
        n.var2 <- length(var2)
        var1 <- rep(var1, times = n.var2)
        res <- sapply(1:n.var2, function(i){
            stats::as.formula(paste(var1[i], var2[i], sep = sep))
        })
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(repVar1 && !missing(var1)){var1 <- rep(var1, length(var2))}
    res <- list(var1 = var1,
                var2 = if(!missing(var2)){var2}else{NULL} 
    )
  }
 
  ## export 
  return(res)
}

## * initVarLinks
#' @rdname initVarLink
#' @export
initVarLinks <- function(var1, format = "list",...){
        
    if("formula" %in% class(var1)){
        res <- initVarLink(var1, repVar1 = TRUE, format = format,
                           ...)
    }else {
        res <- sapply(var1, function(x){
            initVarLink(x, repVar1 = TRUE, format = format,
                        ...)
        })
        if(format == "list"){
            res <- list(var1 = unname(unlist(res["var1",])),
                        var2 = unname(unlist(res["var2",])))
        }else{
            res <- unname(unlist(res))
        }
    }

    return(res)
    
}

## * Documentation - selectResponse
#' @title Response Variable of a Formula
#' @description Return the reponse variable contained in the formula.
#' @name selectResponse
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to lower levels functions.
#'
#' @examples
#' selectResponse(Y1~X1+X2)
#' selectResponse(Y1~X1+X2, type = "vars")
#' selectResponse(Surv(event,time)~X1+X2, type = "vars")
#' 
#' selectResponse(Y1~X1+Y1)
#' selectResponse(Y1+Y2~X1+Y1, type = "vars")
#' 
#' selectResponse(~X1+X2)
#' selectResponse(~X1+X2, type = "vars")
#' 
#' @rdname selectResponse
#' @export
`selectResponse` <-  function(x, ...) UseMethod("selectResponse")

## * selectResponse.formula
#' @rdname selectResponse
#' @method selectResponse formula
#' @export
selectResponse.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[2]]
    if(type == "vars"){
      res <- all.vars(res)
    }
  }else{
    res <- NULL
  }
  
  return(res)
}

## * Documentation - selectRegressor
#' @title Regressor of a Formula.
#' @description Return the regressor variables contained in the formula
#' @name selectRegressor
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to lower levels functions.
#' 
#' @examples
#' 
#' selectRegressor(Y1~X1+X2)
#' selectRegressor(Y1~X1+X2, type = "vars")
#' 
#' selectRegressor(Y1~X1+Y1)
#' selectRegressor(Y1+Y2~X1+Y1, type = "vars")
#' 
#' selectRegressor(~X1+X2)
#' selectRegressor(~X1+X2, type = "vars")
#' 
#' @rdname selectRegressor
#' @export
`selectRegressor` <-  function(x, ...) UseMethod("selectRegressor")

## * selectRegressor.formula
#' @rdname selectRegressor
#' @method selectRegressor formula
#' @export
selectRegressor.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[3]]
    
  }else if(length(x)==2){
    res <- x[[2]]
  }else{
    res <- NULL
  }
  if(type == "vars"){
    res <- all.vars(res)
  }
  
  return(res)
}



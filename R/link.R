## * findNewLink
## ** doc findNewLink
#' @title Find all New Links Between Variables
#' @description Find all new links between variables (copied from lava::modelsearch).
#' 
#' @name findNewLink
#' 
#' @param x a lvm model
#' @param data an optional dataset used to identify the categorical variables if not specified in the lvm object.
#' @param exclude.var all links related to these variables will be ignore.
#' @param rm.latent_latent ignore links relating two latent variables.
#' @param rm.endo_endo ignore links relating two endogenous variables
#' @param rm.latent_endo ignore links relating one endogenous variable and one latent variable
#' @param output return the names of the variables to link ("names") or their position ("index")
#' @param ... additional arguments to be passed to lower levels functions.
#'
#' @return A list
#' 
#' @examples
#' library(lava)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' categorical(m,labels=c("M","F","MF")) <- ~X1
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' findNewLink(m, exclude.var = "X1")
#' 
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' findNewLink(m, rm.endo = TRUE, output = "index")
#' @export
`findNewLink` <-
  function(x, ...) UseMethod("findNewLink")

## ** method findNewLink.lvm
#' @export
#' @rdname findNewLink
findNewLink.lvm <- function(x, data = NULL,
                            exclude.var = NULL, rm.latent_latent= FALSE, rm.endo_endo= FALSE, rm.latent_endo= FALSE,
                            output = "names", ...){

    match.arg(output, choices = c("names","index"))
    if(is.null(data)){        
        data <- lava::sim(x, n = 1)
    }
   
    ## *** convertion to dummy variable name for categorical variables
    xF <- lava_categorical2dummy(x, data)
    AP <- with(lava::index(xF$x), A + t(A) + P)


    if(!is.null(exclude.var)){
        exclude.var <- var2dummy(xF, exclude.var)
    }
    
    if( any(exclude.var %in% colnames(AP) == FALSE) ){
        wrong.var <- exclude.var[exclude.var %in% colnames(AP) == FALSE]
        stop("unknown variable to exclude \n",
             "variable(s): \"",paste(wrong.var, collapse = "\" \""),"\"\n")
        }

    ## *** loop over links
    restricted <- c()
    directional <- c()
    for (i in seq_len(ncol(AP) - 1)){
        for (j in seq(i + 1, nrow(AP))){

            var.i <- rownames(AP)[i]
            var.j <- rownames(AP)[j]
          
            if(!is.null(exclude.var) && (var.i %in% exclude.var || var.j %in% exclude.var)){
                next
            }

            isLatent.i <- var.i %in% latent(xF$x)
            isLatent.j <- var.j %in% latent(xF$x)
            isEndogenous.i <- var.i %in% endogenous(xF$x)
            isEndogenous.j <- var.j %in% endogenous(xF$x)
            isExogenous.i <- var.i %in% exogenous(xF$x)
            isExogenous.j <- var.j %in% exogenous(xF$x)

            if(rm.latent_latent && isLatent.i && isLatent.j){
                next
            }
            if(rm.endo_endo && isEndogenous.i && isEndogenous.j){
                next
            }
            if(isExogenous.i && isExogenous.j){
                next
            }
            if(rm.latent_endo && ( (isLatent.i && isEndogenous.j) || (isEndogenous.i && isLatent.j) )){
                next
            }
      
            if (AP[j, i] == 0){
                restricted <- rbind(restricted, c(i, j))
                directional <- c(directional, (isExogenous.i+isExogenous.j)>0 )
            }
        }
    }

    ## *** export  
    out <- list(M.links = restricted,
                links = NULL,
                directional = directional)

    if(!is.null(restricted)){
        M.names <- cbind(rownames(AP)[restricted[,1]],
                         colnames(AP)[restricted[,2]])
        out$links <- paste0(M.names[,1], lava.options()$symbols[2-directional],M.names[,2])
        if(output == "names"){
            out$M.links <- M.names
        }
    }
   
    return(out)  
}

## * addLink
## ** doc addLink
#' @title Add a New Link Between Two Variables in a LVM
#' @rdname addLink
#' @description Generic interface to add links to lvm objects.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link to be added to the lvm
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param allVars all the existing variables.
#' @param covariance does the link is bidirectional. Ignored if one of the variable is exogenous.
#' @param warnings should a warning be displayed when no link is added
#' @param ... additional arguments to be passed to lower levels functions.
#'
#' @details
#' The argument allVars is useful for \code{lvm.reduce} object where the command \code{vars(x)} does not return all variables. The command \code{vars(x, xlp = TRUE)} must be used instead.
#' 
#' 
#' @examples
#' library(lava)
#' set.seed(10)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' m2 <- m
#' 
#' addLink(m, x1 ~ y1, covariance = FALSE)
#' addLink(m, y1 ~ x1, covariance = FALSE)
#' coef(addLink(m, y1 ~ y2, covariance = TRUE))
#' 
#' addLink(m2, "x1", "y1", covariance = FALSE)
#' addLink(m2, "y1", "x1", covariance = FALSE)
#' newM <- addLink(m, "y1", "y2", covariance = TRUE)
#' coef(newM)
#' @export
`addLink` <-
    function(x, ...) UseMethod("addLink")

## ** method addLink.lvm
#' @export
#' @rdname addLink
addLink.lvm <- function(x,
                        var1,
                        var2,
                        covariance,
                        allVars = vars(x),
                        warnings = FALSE,
                        ...){

    res <- initVarLink(var1, var2, format = "list")
    var1 <- res$var1
    var2 <- res$var2
    
    if(var1 %in% allVars == FALSE){
        if(warnings){
            warning("addLink.lvm: var1 does not match any variable in x, no link is added \n",
                    "var1: ",var1,"\n")
        }
    }
    
    ####
    if(is.na(var2)){
        
        intercept(x) <- stats::as.formula(paste0("~", var1))
        
    }else{
        
        if(var1 == var2){
            if(warnings){
                warning("addLink.lvm: var1 equals var2, no link is added \n",
                        "var1/2: ",var1,"\n")
            }
        }
        
        
        ## if(var2 %in% allVars == FALSE){
        ##     if(warnings){
        ##         warning("addLink.lvm: var2 does not match any variable in x, no link is added \n",
        ##                 "var2: ",var2,"\n")
        ##     }
        ## }
        if(var1 %in% endogenous(x) && var2 %in% endogenous(x)){
            if(missing(covariance)){
                covariance <- TRUE
            }else if(covariance == FALSE){
                stop("addLink.lvm: cannot add a link between two endogenous variable when argument \'covariance\' is FALSE \n")
            }                
        }
    
       
        if(covariance){
            covariance(x) <- stats::as.formula(paste(var1, var2, sep = "~"))  
        }else if(var1 %in% endogenous(x) || var2 %in% exogenous(x)){
            regression(x) <- stats::as.formula(paste(var1, var2,  sep = "~"))
        }else if(var2 %in% endogenous(x) || var1 %in% exogenous(x)){
            regression(x) <- stats::as.formula(paste(var2, var1, sep = "~"))
        }else {
            if(var1 %in% latent(x)){
                regression(x) <- stats::as.formula(paste(var1, var2, sep = "~"))  
            }else if(var2 %in% latent(x)){
                regression(x) <- stats::as.formula(paste(var2, var1, sep = "~"))  
            }else{
                stop("unknow configuration \n")
            }
            
        }
    }
    
    return(x)
}

## ** method addLink.lvm.reduced
#' @rdname addLink
addLink.lvm.reduced <- function(x, ...){
  return(addLink.lvm(x, allVars = vars(x, lp = FALSE, xlp = TRUE) , ...))
}

## * setLink
## ** Documentation - setLink
#' @title Set a Link to a Value
#' @name setLink
#' @description Generic interface to set a value to a link in a lvm object.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param value the value to which the link should be set
#' @param warnings should a warning be displayed when the link is not found in the lvm.
#' @param ... additional arguments to be passed to lower levels functions.
#' 
#' @examples
#' library(lava)
#' set.seed(10)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' m1 <- setLink(m, y3 ~ u, value = 1)
#' estimate(m1, lava::sim(m,1e2))
#' # m1 <- setLink(m, u ~ y3, value = 1)
#' 
#' m2 <- setLink(m, y1 ~ y2, value = 0.5)
#' estimate(m2, lava::sim(m,1e2))
#' @export
`setLink` <-
  function(x, ...) UseMethod("setLink")

## ** method setLink.lvm
#' @rdname setLink
#' @export
setLink.lvm <- function(x, var1, var2, value, warnings = FALSE, ...){

  res <- initVarLink(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### set the link
  if(is.na(var2)){
    intercept(x, stats::as.formula(paste0("~",var1))) <- value
  }else if(paste(var1, var2, sep = "~") %in% stats::coef(x)){
    regression(x, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var1,var2, sep = ",") %in% stats::coef(x)){
    covariance(x, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var2,var1, sep = ",") %in% stats::coef(x)){
    covariance(x, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else{
    if(warnings){
      warning("setLink.lvm: no link was found from var1 to var2, no link is set \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  return(x)
}





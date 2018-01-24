### evalInParentEnv.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (10:48) 
## Version: 
## last-updated: jan 15 2018 (11:45) 
##           By: Brice Ozenne
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Find Object in the Parent Environments
#' 
#' @description Search an object in the parent environments. For internal use.
#' 
#' @param name character string containing the name of the object to get.
#' @param envir the environment from which to look for the object.
#'
#' @keywords internal
evalInParentEnv <- function(name, envir){
  
  ### ** find parent environments
  frames <- sys.status()
  all.frames <- sapply(1:length(frames$sys.frames), function(x){identical(parent.frame(x),globalenv())})
  index.parents <- which(all.frames==FALSE)
  n.parents <- length(index.parents)
  
  ### ** look in parent environments
  iParent <- 1
  res <- NULL
  cv <- FALSE
  while((iParent <= n.parents) && (cv == FALSE)){ # iParent <- 1
    
    res <- try(eval(name, envir = parent.frame(iParent)), silent = TRUE)
      
      if("try-error" %in% class(res)){
        iParent <- iParent + 1     
      }else{
        cv <- TRUE
      }
    
  }
  
  return(res)
}


#----------------------------------------------------------------------
### evalInParentEnv.R ends here

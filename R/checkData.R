### checkData.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 26 2017 (14:25) 
## Version: 
## last-updated: jan 18 2018 (13:29) 
##           By: Brice Ozenne
##     Update #: 23
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - checkData
#' @title Check that Validity of the Dataset
#' @description Check the validity of the dataset used to estimate a lvm.
#' 
#' @name checkData
#' 
#' @param x a lvm model
#' @param data the dataset containing the variables used to estimate the lvm.
#' @param trace [logical] should a message be output to indicate the outcome of the check? 
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#'
#' d <- sim(m,1e2)
#'
#' try(checkData(m, d)) # return an error
#' 
#' checkData(m, d[,-4])
#' 
#' try(checkData(m, d[,-(3:4)])) # return an error
#'
#' @export
`checkData` <-
  function(x, data, trace) UseMethod("checkData")

## * checkData.lvm
#' @rdname checkData
#' @export
checkData.lvm <- function(x, data, trace = TRUE){

    ## ** normalize arguments
    data <- as.data.frame(data)
        
    ## ** check missing names
    vars <- vars(x)
    latent <- latent(x)    
    missingVars <- vars[vars %in% names(data) == FALSE]

    if(length(latent) == 0){
        
        if(length(missingVars)>0){
            if(trace){
                cat("Missing variable in data: ",paste(missingVars, collapse = " "),"\n", sep ="")
            }
            return(invisible(FALSE))
        }else{
            if(trace){
                cat("No issue detected \n")
            }
            return(invisible(TRUE))
        }
        
    }else{
        
        if(!identical(sort(latent),sort(missingVars))){
            if(trace){
                cat("Wrong specification of the latent variables \n",
                    "latent variables according to the LVM: ",paste(latent, collapse = " "),"\n",
                    "missing variables in data: ",paste(missingVars, collapse = " "),"\n", sep = "")
            }
            return(invisible(FALSE))
        }else{
            if(trace){
                cat("No issue detected \n")
            }
            return(invisible(TRUE))
        }
        
    }
}


#----------------------------------------------------------------------
### checkData.R ends here

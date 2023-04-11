## * Documentation
#' @title Extract Data From a Latent Variable Model
#' @description Extract data from a latent variable model.
#' @name extractData
#' 
#' @param object the fitted model.
#' @param design.matrix [logical] should the data be extracted after transformation (e.g. conversion of categorical variables to dummy variables)?
#' Otherwise the original data will be returned.
#' @param as.data.frame [logical] should the output be converted into a \code{data.frame} object?
#' @param rm.na [logical] should the lines containing missing values in the dataset be removed?
#' @param envir [environment] the environment from which to search the data.
#'
#' @return a dataset.
#' @concept extractor
#' @export
`extractData` <-
    function(object, design.matrix, as.data.frame, envir, rm.na){
        UseMethod("extractData", object)
    }

## * Example
#' @rdname extractData
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' n <- 101
#'
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' Id <- findInterval(runif(n), seq(0.1,1,0.1))
#' data.df <- rbind(data.frame(Y=Y1,G="1",Id = Id),
#'            data.frame(Y=Y2,G="2",Id = Id)       
#'            )
#'
#' #### latent variable model ####
#' library(lava)
#' e.lvm <- estimate(lvm(Y ~ G), data = data.df)
#' extractData(e.lvm)
#' extractData(e.lvm, design.matrix = TRUE)
#' 


## * extractData.lvmfit
#' @export
extractData.lvmfit <- function(object, design.matrix = FALSE, as.data.frame = TRUE,
                               envir = environment(), rm.na = TRUE){
    ## ** check arguments
    if(!is.logical(design.matrix)){
        stop("Argument \'design.matrix\' must be of type logical")
    }
    if(!is.logical(as.data.frame)){
        stop("Argument \'as.data.frame\' must be of type logical")
    }

    ## ** extract data
    data <- object$data$model.frame
    if(!inherits(data, "data.frame")){
        data <- as.data.frame(data)
    }

    if(design.matrix){
        keep.cols <- intersect(c("(Intercept)",lava::vars(object)), names(data))
        data <- data[,keep.cols,drop=FALSE]
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** remove missing values relative to the exogenous variables
    test.na <- rowSums(is.na(data[,lava::exogenous(object),drop=FALSE]))
    if(rm.na == TRUE && any(test.na>0)){ ## remove rows corresponding to missing values
        if(!inherits(object,"lvm.missing")){
            data <- data[setdiff(1:NROW(data),which(test.na>0)),,drop=FALSE]
        }else{
            test.na <- rowSums(is.na(data[,exogenous(object)])) > 0
            data <- data[setdiff(1:NROW(data),which(test.na>0)),,drop=FALSE]
            warnings("Missing values in the exogenous variables \n",
                     "May not extract the appropriate dataset \n")
        }
    }

    ## ** export
    return(data)
    
}


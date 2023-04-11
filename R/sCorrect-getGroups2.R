### sCorrect-getGroups2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:58) 
## Version: 
## Last-Updated: Jan 17 2022 (18:43) 
##           By: Brice Ozenne
##     Update #: 170
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Reconstruct the Cluster variable
#' @description Reconstruct the cluster variable.
#' Similar to \code{nlme::getGroups}.
#' @noRd
#'
#' @param object a \code{lvmfit} object.
#' @param data dataset.
#' @param index.Omega  [list] for each cluster, the position of the observed endogenous variables (i.e. how to subset the residual variance-covariance matrix).
#' @param endogenous [character vector] name of the endogenous variables.
#' @param ... [internal] Only used by the generic method.
#'  
#' @return A list containing:
#' \itemize{
#' \item index.cluster: the cluster index for each observation.
#' \item name.cluster: a unique identifier for each cluster.
#' \item n.cluster: the number of clusters.
#' }
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### latent variable model ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' .getGroups2(e.lvm, data = dW)
#' 
#' @concept extractor
#' @keywords internal
`.getGroups2` <-
    function(object, data, index.Omega, endogenous) UseMethod(".getGroups2")

## * .getGroups2.lvm
.getGroups2.lvm <- function(object, data = NULL, index.Omega = NULL, endogenous = NULL){
    if(is.null(data)){
        data <- extractData(object)
    }
    if(is.null(index.Omega)){
        index.Omega <- .getIndexOmega(object, data = data)
    }
    if(is.null(endogenous)){
        endogenous <- lava::endogenous(object)
    }
    n.endogenous <- length(endogenous)

    ## ** find clusters
    n.cluster <- NROW(data)
    name.cluster <- 1:n.cluster
    missing <- any(is.na(index.Omega))
    index.cluster <- unlist(lapply(name.cluster, rep, times = n.endogenous))

    index.Omega <- tapply(index.Omega, index.cluster, function(iVec){list(stats::na.omit(iVec))})
    Uindex.Omega <- unique(index.Omega[sapply(index.Omega,length)>0])
    return(list(index.cluster = index.cluster,
                name.cluster = name.cluster,
                n.cluster = n.cluster,
                index.Omega = index.Omega,
                index2endogenous = stats::setNames(as.list(Uindex.Omega),Uindex.Omega)
                ))
    
}

## * .getGroups2.lvmfit
.getGroups2.lvmfit <- .getGroups2.lvm








######################################################################
### sCorrect-getGroups2.R ends here

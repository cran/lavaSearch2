### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: Jan 17 2022 (23:21) 
##           By: Brice Ozenne
##     Update #: 2407
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - score2
#' @title  Score With Small Sample Correction
#' @description  Extract the (individual) score a the latent variable model.
#' Similar to \code{lava::score} but with small sample correction.
#' @name score2
#'
#' @param object,x a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param indiv [logical] If \code{TRUE}, the score relative to each observation is returned. Otherwise the total score is returned.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the confidence intervals.
#'
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#'
#' @return When argument indiv is \code{TRUE}, a matrix containing the score relative to each sample (in rows)
#' and each model coefficient (in columns). Otherwise a numeric vector of length the number of model coefficients.
#' 
#' @examples
#' #### simulate data ####
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- Sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' #### linear models ####
#' e.lm <- lm(Y~X1+X2+X3, data = d)
#' 
#' #### latent variable models ####
#' m.lvm <- lvm(formula.lvm)
#' e.lvm <- estimate(m.lvm,data=d)
#' e2.lvm <- estimate2(m.lvm,data=d)
#' score.tempo <- score(e2.lvm, indiv = TRUE)
#' colSums(score.tempo)
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`score2` <-
  function(object, indiv, cluster, as.lava, ...) UseMethod("score2")

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit <- function(object, indiv = FALSE, cluster = NULL, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(lava::score(estimate2(object, ssc = ssc, ...), indiv = indiv, cluster = cluster, as.lava = as.lava))

}

## * score2.lvmfit2
#' @rdname score2
#' @export
score2.lvmfit2 <- function(object, indiv = FALSE, cluster = NULL, as.lava = TRUE, ...){
    
    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }
    ## ** define cluster
    if(length(cluster) == 1 && (is.numeric(cluster) || is.character(cluster) || is.factor(cluster))){
        data <- object$sCorrect$data
        if(length(cluster)==1){                
            if(cluster %in% names(data) == FALSE){
                stop("Invalid \'cluster\' argument \n",
                     "Could not find variable \"",cluster,"\" in argument \'data\' \n")
            }
            cluster <- data[[cluster]]
        }
        cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))            
        n.cluster <- length(unique(cluster.index))
    }else if(is.vector(cluster)){
        cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))
        n.cluster <- length(unique(cluster.index))
    }else if(is.null(cluster)){ ## NOTE: cluster is a function in the survival package
        cluster <- NULL
        n.cluster <- object$sCorrect$cluster$n.cluster
        cluster.index <- 1:n.cluster
    }else{
        stop("Do not know how to handle argument cluster of class ",class(cluster),"\n")
    }

    ## ** get score
    score <- object$sCorrect$score
    if(!is.null(cluster)){ ## aggregate score by cluster
        score <- rowsum(score, group = cluster.index, reorder = FALSE)
    }
    
    ## ** export
    score <- score[,names(object$sCorrect$skeleton$originalLink2param),drop=FALSE]
    if(as.lava==FALSE){
        colnames(score) <- as.character(object$sCorrect$skeleton$originalLink2param)
    }
    if(!is.null(cluster)){
        index2.cluster <- tapply(1:length(cluster),cluster,list)
        attr(score,"cluster") <- names(index2.cluster)
    }
    
    if(indiv){
        return(score)
    }else{
        return(colSums(score))
    }
}

## * score.lvmfit2
#' @rdname score2
#' @export
score.lvmfit2 <- function(x, indiv = FALSE, cluster = NULL, as.lava = TRUE, ...){## necessary as first argument of score must be x 
    score2(x, indiv = indiv, cluster = cluster, as.lava = as.lava, ...)
}

## * .score2
#' @title Compute the Corrected Score.
#' @description Compute the corrected score.
#' @name score2-internal
#' 
#' @param n.cluster [integer >0] the number of observations.
#' 
#' @keywords internal
.score2 <- function(dmu, dOmega, epsilon, OmegaM1,
                    missing.pattern, unique.pattern, name.pattern,
                    name.param, name.meanparam, name.varparam,
                    n.cluster, weights){
    if(lava.options()$debug){cat(".score2\n")}

    ## ** Prepare
    out.score <- matrix(NA, nrow = n.cluster, ncol = length(name.param),
                        dimnames = list(NULL,name.param))
    n.pattern <- length(name.pattern)
    
    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iEpsilon.OmegaM1 <- epsilon[iIndex,iY,drop=FALSE] %*% iOmegaM1
        out.score[iIndex,] <- 0 ## initialize (keep NA for missing values)

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- "Y3~eta"
            out.score[iIndex,iP] <- out.score[iIndex,iP] + rowSums(dmu[[iP]][iIndex,iY,drop=FALSE] * iEpsilon.OmegaM1)
        }
        
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.varparam){ # iP <- "eta~~eta"
            term2 <- - 1/2 * tr(iOmegaM1 %*% dOmega[[iP]][iY,iY,drop=FALSE])            
            term3 <- 1/2 * rowSums(iEpsilon.OmegaM1 %*% dOmega[[iP]][iY,iY,drop=FALSE] * iEpsilon.OmegaM1)
            out.score[iIndex,iP] <- out.score[iIndex,iP] + as.double(term2) + term3
        }        
    }

    ## ** export
    if(!is.null(weights)){
        out.score <- sweep(out.score, STATS = weights, MARGIN = 1, FUN = "*")
    }
    return(out.score)
}


#----------------------------------------------------------------------
### score2.R ends her

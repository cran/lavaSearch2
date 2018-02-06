### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: feb  5 2018 (17:26) 
##           By: Brice Ozenne
##     Update #: 2201
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - score2
#' @title Compute the Score Directly from the Gaussian Density
#' @description Compute the score directly from the gaussian density.
#' @name score2
#' 
#' @param object a linear model or a latent variable model.
#' @param p [numeric vector, optional] vector of coefficients at which to evaluate the score.
#' @param data [data.frame, optional] data set.
#' @param cluster [vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param indiv [logical] should the score relative to each observation be exported?
#' Otherwise the total score (i.e. sum over all observations) will be exported.
#' @param return.vcov.param [logical] should the variance covariance matrix of the coefficients be included in the output?
#' @param ... arguments to be passed to \code{\link{residuals2}}
#'
#' @return A matrix.
#' 
#' @details The log-likelihood of a lvm can written:
#' \deqn{
#'   l(\theta|Y,X) = \sum_{i=1}^{n} - \frac{p}{2} log(2\pi) - \frac{1}{2} log|\Omega(\theta))| - \frac{1}{2} (Y_i-\mu_i(\theta)) \Omega^{-1} (Y_i-\mu_i(\theta))^t
#' }
#' So the score is:
#' \deqn{
#'   s(\theta|Y,X) = \sum_{i=1}^{n}
#' - \frac{1}{2} tr(\Omega(\theta)^{-1} \frac{\partial \Omega(\theta)}{\partial \theta})
#' +  \frac{\partial \mu_i(\theta)}{\partial \theta} \Omega^{-1} (Y_i-\mu_i(\theta))^t 
#' + \frac{1}{2} (Y_i-\mu_i(\theta)) \Omega^{-1} \frac{\partial \Omega(\theta)}{\partial \theta} \Omega^{-1} (Y_i-\mu_i(\theta))^t
#' }
#' with:
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \nu} = 1 }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial K} = X_i }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \alpha} = (1-B)^{-1}\Lambda }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \Gamma} = X_i(1-B)^{-1}\Lambda }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \lambda} = (\alpha + X_i \Gamma)(1-B)^{-1}\delta_{\lambda \in \Lambda} }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial b} = (\alpha + X_i \Gamma)(1-B)^{-1}\delta_{b \in B}(1-B)^{-1}\Lambda }
#' and:
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \psi} = \Lambda^t (1-B)^{-t} \delta_{\psi \in \Psi} (1-B) \Lambda }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \sigma} = \delta_{\sigma \in \Sigma} }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \lambda} = \delta_{\lambda \in \Lambda} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda + \Lambda^t (1-B)^{-t} \Psi (1-B)^{-1} \delta_{\lambda \in \Lambda} }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial b} = \Lambda^t (1-B)^{-t} \delta_{b \in B} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda - \Lambda^t (1-B)^{-t} \Psi (1-B)^{-1} \delta_{b \in B} (1-B)^{-1} \Lambda}
#'
#' @examples
#'
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m,sim(m,1e2))
#' score2(e)
#'
#' @concept small sample inference
#' @export
`score2` <-
  function(object, ...) UseMethod("score2")


## * score2.lm
#' @rdname score2
#' @export
score2.lm <- function(object,
                      indiv = TRUE, return.vcov.param = FALSE,
                      ...){
  
    epsilon <- residuals2(object, return.vcov.param = TRUE, ...)
    vcov.param <- attr(epsilon, "vcov.param")
    sigma2 <- attr(epsilon, "sigma2") 
    attr(epsilon, "vcov.param") <- NULL
    attr(epsilon, "sigma2") <- NULL

    ## ** compute score
    X <- stats::model.matrix(object)
    out.score <- cbind(sweep(X, MARGIN = 1, FUN = "*",  STATS = epsilon) / sigma2,
                       sigma2 = -1/(2*sigma2) + 1/(2*sigma2^2) * epsilon^2)

    if(indiv == FALSE){
        out.score <- colSums(out.score)
    }

    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param
    }
    return(out.score)

}


## * score2.gls
#' @rdname score2
#' @export
score2.gls <- function(object, cluster, p = NULL, data = NULL,
                       indiv = TRUE, return.vcov.param = FALSE,
                       ...){

### ** Compute residuals and partial derivatives
    epsilon <- residuals2(object, cluster = cluster, p = p, data = data,
                          return.vcov.param = return.vcov.param,
                          return.prepareScore2 = TRUE, ...)

    OPS2 <- attr(epsilon, "prepareScore2")
    vcov.param <- attr(epsilon, "vcov.param")
    attr(epsilon, "prepareScore2") <- NULL
    attr(epsilon, "vcov.param") <- NULL

### ** Compute score
    out.score <- .score2(dmu.dtheta = OPS2$dmu.dtheta,
                         dOmega.dtheta = OPS2$dOmega.dtheta,
                         epsilon = epsilon,
                         Omega = OPS2$Omega,                         
                         ls.indexOmega = OPS2$ls.indexOmega,
                         indiv = indiv,
                         name.param = OPS2$name.param,
                         n.param = OPS2$n.param,
                         n.cluster = OPS2$n.cluster)     
    
### ** Export
    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param
    }
    return(out.score)
}

## * score2.lme
#' @rdname score2
#' @export
score2.lme <- score2.gls

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit <- function(object, p = NULL, data = NULL, 
                          indiv = TRUE, return.vcov.param = FALSE,
                          ...){

### ** Compute residuals
    epsilon <- residuals2(object,
                          p = p,
                          data = data,
                          return.vcov.param = return.vcov.param,
                          return.prepareScore2 = TRUE, ...)

    OPS2 <- attr(epsilon, "prepareScore2")
    vcov.param <- attr(epsilon, "vcov.param")
    attr(epsilon, "prepareScore2") <- NULL
    attr(epsilon, "vcov.param") <- NULL

### ** Compute score
    out.score <- .score2(dmu.dtheta = OPS2$dtheta$dmu.dtheta,
                         dOmega.dtheta = OPS2$dtheta$dOmega.dtheta,
                         epsilon = epsilon,
                         Omega = OPS2$Omega,
                         ls.indexOmega = NULL,
                         indiv = indiv,                         
                         name.param = OPS2$name.param,
                         n.param = OPS2$n.param,
                         n.cluster = OPS2$n.cluster)     
    
### ** Export
    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param[OPS2$name.param, OPS2$name.param, drop=FALSE]
    }
    return(out.score)
}




## * .score2
.score2 <- function(dmu.dtheta, dOmega.dtheta, epsilon,
                    Omega, ls.indexOmega,
                    indiv, 
                    name.param, n.param, n.cluster){

    clusterSpecific <- !is.null(ls.indexOmega)
    name.meanparam <- names(dmu.dtheta)
    name.vcovparam <- names(dOmega.dtheta)
    out.score <- matrix(0, nrow = n.cluster, ncol = n.param,
                        dimnames = list(NULL,name.param))

### ** Individual specific Omega (e.g. presence of missing values)
    if(clusterSpecific){
        for(iC in 1:n.cluster){
            iOmega.tempo <- chol2inv(chol(Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]))
            epsilon.iOmega.tempo <- iOmega.tempo %*% cbind(epsilon[iC,ls.indexOmega[[iC]]])

            ## *** Compute score relative to the mean coefficients
            for(iP in name.meanparam){ # iP <- name.meanparam[1]
                out.score[iC,iP] <- out.score[iC,iP] + dmu.dtheta[[iP]][iC,ls.indexOmega[[iC]]] %*% epsilon.iOmega.tempo
            }

            ## *** Compute score relative to the variance-covariance coefficients
            for(iP in name.vcovparam){ # iP <- name.vcovparam[1]
                dOmega.dtheta.tempo <-  dOmega.dtheta[[iP]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                                                            
                term2 <- - 1/2 * tr(iOmega.tempo %*% dOmega.dtheta.tempo)
                term3 <- 1/2 * sum(epsilon.iOmega.tempo * dOmega.dtheta.tempo %*% epsilon.iOmega.tempo)
                out.score[iC,iP] <- out.score[iC,iP] + as.double(term2) + term3 
            }
        }
    }
            
    ### ** Same for all individuals
    if(clusterSpecific == FALSE){
        iOmega <- chol2inv(chol(Omega))
        epsilon.iOmega <- epsilon %*% iOmega

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- 1
            out.score[,iP] <- out.score[,iP] + rowSums(dmu.dtheta[[iP]] * epsilon.iOmega)            
        }
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.vcovparam){ # iP <- 1
            term2 <- - 1/2 * tr(iOmega %*% dOmega.dtheta[[iP]])            
            term3 <- 1/2 * rowSums(epsilon.iOmega %*% dOmega.dtheta[[iP]] * epsilon.iOmega)
            out.score[,iP] <- out.score[,iP] + as.double(term2) + term3 
        }        
    }

### ** export
    if(indiv==FALSE){
        out.score <- colSums(out.score)
    }
    return(out.score)
}

#----------------------------------------------------------------------
### score2.R ends her

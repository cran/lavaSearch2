### matrixPower.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 23 2017 (16:52) 
## Version: 
## last-updated: feb  5 2018 (16:45) 
##           By: Brice Ozenne
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

##' @title Power of a Matrix
##' @description Compute the power of a matrix.
##' 
##' @param object a matrix.
##' @param power [numeric] power to be applied to the matrix.
##' @param symmetric [logical] is the matrix symmetric? Argument passed to the function \code{eigen}.
##' @param tol [numeric >0] the threshold under which the eigenvalues are set to 0.
##' 
##' @return A matrix.
##' 
##' @examples
##' matrixPower <- lavaSearch2:::matrixPower
##' 
##' M <- matrix(rnorm(2e2),20,10)
##' Sigma <- var(M)
##' Sigma.half <- matrixPower(Sigma, power = 1/2, symmetric = FALSE)
##' round(crossprod(Sigma.half)-Sigma,5)
##' 
##' Sigma.m1 <- matrixPower(Sigma, power = -1, symmetric = FALSE)
##' round(Sigma.m1 %*% Sigma,5)
##' 
##' @keywords internal
matrixPower <- function(object, power, symmetric, tol = 1e-12){
    object.eigen <- eigen(object, symmetric = symmetric)
    object.eigen$values[abs(object.eigen$values) < tol] <- 0 ## abs is not really necessary since the eigenvalues are positives
    object.eigen$values <- object.eigen$values^power
    
    return(object.eigen$vectors %*% sweep(t(object.eigen$vectors),
                                     MARGIN = 1,
                                     FUN = "*",
                                     STATS = object.eigen$values)
           )
}

#----------------------------------------------------------------------
### matrixPower.R ends here

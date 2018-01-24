### matrixPower.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 23 2017 (16:52) 
## Version: 
## last-updated: jan 17 2018 (19:05) 
##           By: Brice Ozenne
##     Update #: 27
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
##' @param x a matrix.
##' @param power power to be applied to the matrix.
##' @param symmetric is the matrix symmetric? Argument passed to the function \code{eigen}.
##' @param tol the threshold under which the eigenvalues are set to 0.
##' @return a matrix.
##'
##' @examples
##' M <- matrix(rnorm(2e2),20,10)
##' Sigma <- var(M)
##' Sigma.half <- matrixPower(Sigma, power = 1/2, symmetric = FALSE)
##' round(crossprod(Sigma.half)-Sigma,5)
##' 
##' Sigma.m1 <- matrixPower(Sigma, power = -1, symmetric = FALSE)
##' round(Sigma.m1 %*% Sigma,5)
##' 
##' @export
matrixPower <- function(x, power, symmetric, tol = 1e-12){
    x.eigen <- eigen(x, symmetric = symmetric)
    x.eigen$values[abs(x.eigen$values) < tol] <- 0 ## abs is not really necessary since the eigenvalues are positives
    x.eigen$values <- x.eigen$values^power
    return(x.eigen$vectors %*% sweep(t(x.eigen$vectors),
                                     MARGIN = 1,
                                     FUN = "*",
                                     STATS = x.eigen$values)
           )
}

#----------------------------------------------------------------------
### matrixPower.R ends here

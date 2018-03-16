### matrixPower.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 23 2017 (16:52) 
## Version: 
## last-updated: mar  7 2018 (12:17) 
##           By: Brice Ozenne
##     Update #: 48
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
##' ## symmetric matrix
##' set.seed(10)
##' M <- matrix(rnorm(20*6),20,6)
##' Sigma <- var(M)
##' Sigma.half <- matrixPower(Sigma, power = 1/2, symmetric = TRUE)
##' round(Sigma.half %*% Sigma.half - Sigma,5)
##' 
##' iSigma <- matrixPower(Sigma, power = -1, symmetric = TRUE)
##' round(iSigma %*% Sigma,5)
##' 
##' iSigma.half <- matrixPower(Sigma, power = -1/2, symmetric = TRUE)
##' round(iSigma.half %*% iSigma.half - iSigma,5)
##' 
##' ## non symmetric matrix
##' set.seed(10)
##' M <- matrix(abs(rnorm(9)), 3, 3) + diag(1,3,3)
##' M-t(M)
##' 
##' iM <- matrixPower(M, power = -1, symmetric = FALSE)
##' round(iM %*% M,5)
##' 
##' iM.half <- matrixPower(M, power = -1/2, symmetric = FALSE)
##' round(iM.half %*% iM.half %*% M,5)
##' 
##' @keywords internal
matrixPower <- function(object, power, symmetric, tol = 1e-12){
    object.eigen <- eigen(object, symmetric = symmetric)

    if(power<0 && any(object.eigen$values<tol)){
        warning("negative eigenvalues are set to ",tol,"\n")
        object.eigen$values[object.eigen$values < tol] <- tol
    }
    nRow <- length(object.eigen$values)
    D <- diag(object.eigen$values^power, nrow = nRow, ncol = nRow)

    if(symmetric){
        return(object.eigen$vectors %*% D %*% t(object.eigen$vectors))
    }else{
        return(object.eigen$vectors %*% D %*% solve(object.eigen$vectors))
    }

}

#----------------------------------------------------------------------
### matrixPower.R ends here

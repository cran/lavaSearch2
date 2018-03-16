### test1-sCorrect-adjustedResiduals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:21) 
## Version: 
## Last-Updated: mar 13 2018 (13:25) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
context("sCorrect: small sample correction")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]

## * linear regression [lm,gls,lvm]
## ** model fit and sCorrect
e.lvm <- estimate(lvm(Y1~X1+X2+Gender), data = d)
e.lm <- lm(Y1~X1+X2+Gender, data = d)
e.gls <- gls(Y1~X1+X2+Gender, data = d, method = "ML")

e2.lvm <- e.lvm
e2.gls <- e.gls
e2.lm <- e.lm

sCorrect(e2.lvm) <- TRUE
sCorrect(e2.gls, cluster = 1:n) <- TRUE
sCorrect(e2.lm) <- TRUE

## ** test adjusted residuals
test_that("residuals2 match residuals.lm (lm adjusted)", {
    epsilon.lm <- residuals(e.lm)
    X <- model.matrix(e.lm, d)
    iH <- diag(1,n,n) - X %*% solve(t(X) %*% X) %*% t(X)
    GS1 <- epsilon.lm/diag(iH)^(1/2)

    expect_equal(as.double(e2.lm$sCorrect$epsilon),
                 as.double(GS1))
    expect_equal(as.double(e2.lvm$sCorrect$epsilon),
                 as.double(GS1))
    expect_equal(as.double(e2.gls$sCorrect$epsilon),
                 as.double(GS1))
})



## * multivariate linear models
## ** model fit
ls.lm <- list(lm(Y1~X1,d),lm(Y2~X2,d),lm(Y3~X1+X3,d))
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)

sCorrect(e.lvm) <- TRUE

test_that("residuals2 match residuals.lm", {

    epsilon.lm <- lapply(ls.lm, residuals)
    X <- lapply(ls.lm, model.matrix)
    iH <- lapply(X, function(x){
        diag(1,NROW(x),NROW(x)) - x %*% solve(t(x) %*% x) %*% t(x)
    })
    GS <- mapply(iH,epsilon.lm, FUN = function(x,y){
        y/diag(x)^(1/2)
    })
    expect_equal(as.double(e.lvm$sCorrect$epsilon),
                 as.double(GS))   
})



##----------------------------------------------------------------------
### test1-sCorrect-adjustedResiduals.R ends here

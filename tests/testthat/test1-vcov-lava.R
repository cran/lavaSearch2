### test-vcov.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (11:44) 
## Version: 
## last-updated: jan 19 2018 (14:58) 
##           By: Brice Ozenne
##     Update #: 66
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * header
if(FALSE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

library(pbkrtest)
library(lme4)
lava.options(symbols = c("~","~~"))

context("vcov")
n <- 5e1

## * Model vcov
## ** function
rmAttr <- function(x, name.rm = NULL, name.keep){
    if(is.null(name.rm)){
        name.rm <- names(attributes(x))
    }
    for(iAttr in name.rm){
        attr(x, iAttr) <- NULL
    }
    return(x)
}

## ** linear regression
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("linear regression (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("linear regression (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("linear regression (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("linear regression: constrains",{
    m <- lvm(Y[0:2]~X1+1*X2)
    e <- estimate(m, d)
    
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))  
    
    m <- lvm(Y~beta*X1+beta*X2)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))  
})


## ** multiple linear regression
## *** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("multiple linear regression (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))

    expect_equal(unname(test), unname(GS))    
})

test_that("multiple linear regression (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("multiple linear regression (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("multiple linear regressions: constrains",{
    m <- lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

## *** with covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("multiple linear regression, covariance link (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("multiple linear regression, covariance link (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("multiple linear regression, covariance link (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

## ** latent variable model
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## *** factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("factor model (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("factor model (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("factor model (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("factor model: fixed coefficients",{
    m <- lvm(Y1~1*eta+1*X2,Y2~1*eta,Y3~1*eta)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("factor model: constrains",{
    m <- lvm(Y1~1*eta+X2,Y2~lambda*eta+X2,Y3~lambda*eta,eta ~ beta*X2+beta*X1)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})


## *** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("2 factor model: constrains",{
    m <- lvm(Y1~1*eta1+X2,Y2~lambda*eta1+X2,Y3~lambda*eta1,eta1 ~ beta*X2+beta*X1,
             Z1~0+eta2,Z2~lambda*eta2,Z3~eta2)
    e <- estimate(m, d)

     test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

## *** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model, covariance (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model, covariance (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model, covariance (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

## *** 2 factor model (correlation LV)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model, correlation (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model, correlation (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model, correlation (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

## * Corrected vcov
n <- 2e1

## ** linear model
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.lvm <- estimate(lvm(Y~X1+X2+X3),d)
e.lm <- lm(Y~X1+X2+X3, data = d)
Sigma.lvm <- attr(residuals2(e.lvm, adjust.residuals = TRUE, return.vcov.param = TRUE), "vcov.param")
Sigma.lm <- attr(residuals2(e.lm, adjust.residuals = TRUE, return.vcov.param = TRUE), "vcov.param")

test_that("Corrected vcov - linear model",{
    expect_equal(unname(Sigma.lvm), unname(Sigma.lm))

    Sigma.uncorrected <- vcov(e.lvm)
    attr(Sigma.uncorrected, "det") <- NULL
    attr(Sigma.uncorrected, "pseudo") <- NULL
    attr(Sigma.uncorrected, "minSV") <- NULL

    p <- 4
    factor_beta <- (n+p)/n 
    expect_equal(Sigma.lvm[1:p,1:p], factor_beta * Sigma.uncorrected[1:p,1:p])

    factor_sigma <- n/(n-p)*(n+p)^2/n^2
    expect_equal(Sigma.lvm[p+1,p+1], factor_sigma * Sigma.uncorrected[p+1,p+1])
})

## ** multiple linear model

## *** simulation
m <- lvm(c(Y1,Y2)~X1+X2)
set.seed(10)
df <- sim(m,2e1)

## *** model fit
e1.lm <- lm(Y1~X1+X2, data = df)
e2.lm <- lm(Y2~X1+X2, data = df)
e.lvm <- estimate(m, data = df)

## *** tests
Sigma.lvm <- attr(residuals2(e.lvm, return.vcov.param = TRUE, adjust.residuals = TRUE),
                 "vcov.param")

Sigma.lvm[c("Y1","Y1~X1","Y1~X2"),c("Y1","Y1~X1","Y1~X2")]/vcov(e1.lm)
Sigma.lvm[c("Y2","Y2~X1","Y2~X2"),c("Y2","Y2~X1","Y2~X2")]/vcov(e2.lm)

Sigma.lvm[c("Y1~~Y1"),c("Y1~~Y1")]/(2*sigma(e1.lm)^4/(NROW(df)-length(coef(e1.lm))))
Sigma.lvm[c("Y2~~Y2"),c("Y2~~Y2")]/(2*sigma(e2.lm)^4/(NROW(df)-length(coef(e2.lm))))



## ** mixed model

## *** simulate
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender,X1~1,X2~1))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,
                     id.vars = c("G","Id","Gender","X1","X2"),
                     variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]

## *** model fit
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender))
e.lvm <- estimate(m, dW)


e.lmer <- lme4::lmer(value ~ time + G + Gender + (1|Id),
                     data = dL, REML = FALSE)


## *** tests
Sigma0.lvm <- vcov(e.lvm)
SigmaAdj.lvm <- attr(residuals2(e.lvm, adjust.residuals = TRUE, return.vcov.param = TRUE), "vcov.param")
## Sigma0.lvm/SigmaAdj.lvm
SigmaAdjRed.lvm <- SigmaAdj.lvm[1:5,1:5]
Sigma.GS <- vcovAdj(e.lmer)

## not far from KR correction
SigmaAdjRed.lvm/Sigma.GS
Sigma0.lvm[1:5,1:5]/Sigma.GS


#----------------------------------------------------------------------
### test-vcov.R ends here





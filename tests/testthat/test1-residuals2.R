### test-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:08) 
## Version: 
## Last-Updated: feb  5 2018 (17:17) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
if(FALSE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

library(nlme)
lava.options(symbols = c("~","~~"))

context("residuals2")

## * Simulation
n <- 5e1

mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]
dL$Z1 <- rnorm(NROW(dL))

## * raw residuals
## ** univariate linear model
m <- lvm(Y~X)
d <- sim(m,1e2)

test_that("residuals2 match residuals.lm (single lm)", {
    GS <- residuals(lm(Y~X, data = d))
    e <- estimate(lvm(Y~X), d)
    res <- residuals(e)
    res2 <- residuals2(e, bias.correct = FALSE)

    expect_equal(as.double(res),as.double(GS))
    expect_equal(as.double(res2),as.double(GS))
})


## ** multivariate linear models
m <- lvm(Y~G+X,G~X)
d <- sim(m,1e2)

test_that("residuals2 match residuals.lm (multiple lm)", {
    GS <- cbind(residuals(lm(Y~1, data = d)),
                residuals(lm(G~1, data = d)))
    e <- estimate(lvm(Y~1,G~1), d)
    res2 <- residuals2(e, bias.correct = FALSE)
    res <- residuals(e)

    expect_equal(unname(res),unname(GS))
    expect_equal(unname(res2),unname(GS))

    GS <- cbind(residuals(lm(Y~G+X, data = d)),
                residuals(lm(G~1, data = d)))
    e <- estimate(lvm(Y~G+X,G~1), d)

    res2 <- residuals2(e, bias.correct = FALSE)
    res <- residuals(e)

    
    ## expect_equal(as.double(res),as.double(GS))
    ## note: vcov(lm(Y~G+X, data = d))/vcov(e)[c("Y","Y~G","Y~X"),c("Y","Y~G","Y~X")]
    
    expect_equal(unname(coef(e)[c("Y","Y~G","Y~X")]),
                 unname(coef(lm(Y~G+X,data=d))), tol = 1e-5)
    expect_equal(res2[,"Y"],unname(GS[,1]), tol = 1e-4)
    expect_equal(res2[,"G"],unname(GS[,2]))
    
    GS <- cbind(residuals(lm(Y~G+X, data = d)),
                residuals(lm(G~X, data = d)))
    e <- estimate(lvm(Y~G+X,G~X), d)

    res2 <- residuals2(e, bias.correct = FALSE)
    res <- residuals(e)

    ## expect_equal(as.double(res),as.double(GS))
    ## note: vcov(lm(Y~G+X, data = d))/vcov(e)[c("Y","Y~G","Y~X"),c("Y","Y~G","Y~X")]
    expect_equal(as.double(res2),as.double(GS))
})

## ** mixed model
## *** versus nlme
mSim <- lvm(c(Y1~1*eta1,Y2~1*eta1,Y3~1*eta1,eta1~G1))
latent(mSim) <- ~eta1
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
dW <- sim(mSim, 5e1, latent = FALSE)
dL <- reshape2::melt(dW, id.vars = c("Id","G1"))


test_that("equivalence residuals2.lvm residuals.lvm", {
    
    m <- lvm(c(Y1[mu1:sigma]~1*eta1,Y2[mu2:sigma]~1*eta1,Y3[mu3:sigma]~1*eta1,eta1~G1))
    latent(m) <- ~eta1
    e.lvm <- estimate(m,dW)

    e.gls <- nlme::gls(value ~ variable + G1, data = dL,
                       correlation = corCompSymm(form =~ variable|Id),
                       method = "ML")
    e.lme <- nlme::lme(value ~ variable + G1, data = dL,
                       random =~ 1|Id,
                       method = "ML")

    expect_equal(as.double(logLik(e.lvm)),as.double(logLik(e.gls)))
    expect_equal(as.double(logLik(e.lvm)),as.double(logLik(e.lme)))
    
    test.gls <- residuals2(e.gls, bias.correct = FALSE)
    test.lme <- residuals2(e.lme, bias.correct = FALSE)
    expect_equal(unname(test.gls),unname(test.lme))
    
    test.lvm <- residuals2(e.lvm, bias.correct = FALSE)
    expect_equal(unname(test.lvm),unname(test.gls))

    GS.lme <- as.double(residuals(e.lme, type = "response", level = 0))
    GS.gls <- as.double(residuals(e.gls))
    expect_equal(GS.lme,GS.gls)
    
    expect_equal(GS.lme,as.double(test.gls))
    expect_equal(GS.lme,as.double(residuals(e.lvm)))
})

## *** versus lvm
m <- lvm(c(Y1~1*eta1,Y2~1*eta1,Y3~1*eta1,eta1~beta*G1,
           Z1~1*eta2,Z2~1*eta2,Z3~1*eta2,eta2~beta*G2)
         )
latent(m) <- ~eta1+eta2
d <- sim(m, 5e1)
e.lvm <- estimate(m,d)

test_that("equivalence residuals2.lvm residuals.lvm", {
    test <- residuals2(e.lvm, bias.correct = FALSE)
    GS <- residuals(e.lvm)
    expect_equal(GS,test)
})

## * adjusted residuals

## ** univariate linear model
m <- lvm(Y~X)
n <- 1e2
d <- sim(m,n)

test_that("residuals2 match residuals.lm (lm adjusted)", {
    e.lm <- lm(Y~X, data = d)
    epsilon.lm <- residuals(e.lm)
    X <- model.matrix(e.lm, d)
    iH <- diag(1,n,n) - X %*% solve(t(X) %*% X) %*% t(X)
    GS1 <- epsilon.lm/diag(iH)^(1/2)
    
    e.lvm <- estimate(lvm(Y~X), d)
    res2 <- residuals2(e.lvm, bias.correct = TRUE)

    expect_equal(as.double(res2),as.double(GS1))
})



## ** multivariate linear models
m <- lvm(Y~G+X,G~X)
n <- 1e2
d <- sim(m,n)

test_that("residuals2 match residuals.lm", {
    ## first model
    e.lm1 <- lm(Y~G+X, data = d)
    res1 <- residuals2(e.lm1, bias.correct = TRUE)

    epsilon.lm1 <- residuals(e.lm1)
    X1 <- model.matrix(e.lm1, d)
    iH1 <- diag(1,n,n) - X1 %*% solve(t(X1) %*% X1) %*% t(X1)

    expect_equal(as.double(epsilon.lm1/diag(iH1)^(1/2)),
                 as.double(res1))

    ## second model
    e.lm2 <- lm(G~X, data = d)
    res2 <- residuals2(e.lm2, bias.correct = TRUE)
    epsilon.lm2 <- residuals(e.lm2)
    X2 <- model.matrix(e.lm2, d)
    iH2 <- diag(1,n,n) - X2 %*% solve(t(X2) %*% X2) %*% t(X2)
    
    expect_equal(as.double(epsilon.lm2/diag(iH2)^(1/2)),
                 as.double(res2))

    ## global
    e.lvm <- estimate(m, d)    
    resTest <- residuals2(e.lvm, bias.correct = TRUE)

    expect_equal(as.double(resTest[,1]),as.double(res1))
    expect_equal(as.double(resTest[,2]),as.double(res2))   
})

##----------------------------------------------------------------------
### test-residuals.R ends here

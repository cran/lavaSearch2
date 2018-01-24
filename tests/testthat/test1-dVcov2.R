### test-dVcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (15:17) 
## Version: 
## Last-Updated: jan 19 2018 (15:55) 
##           By: Brice Ozenne
##     Update #: 65
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
    rm(list = ls(all.names = TRUE))
    library(testthat)
    library(lavaSearch2)
}

library(nlme)
library(lme4)
lava.options(symbols = c("~","~~"))

context("dVcov2")

## * Simulation
n <- 5e1
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender,X1~1,X2~1))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
set.seed(10)
dW <- sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id","Gender","X1","X2"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]

## * linear regression

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2), data = dW)
e.lvm$prepareScore2 <- prepareScore2(e.lvm, second.order = TRUE, usefit = FALSE)
e.gls <- gls(Y1~X1+X2, data = dW, method = "ML")

test_that("linear regression: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, cluster = dW$Id, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, cluster = dW$Id, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.gls, res.gls)
    
})

## * mixed model

## ** Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender)) 
e.lvm <- estimate(m, dW)

e.lmer <- lmer(value ~ time + G + Gender + (1|Id),
               data = dL, REML = FALSE)

e.lme <- lme(value ~ time + G + Gender, random = ~ 1|Id, data = dL, method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corCompSymm(form=~ 1|Id),
             data = dL, method = "ML")

expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed model CS: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.gls, res.gls)

    ## lme
    GS.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lme, res.lme)
})

## ** Unstructured with weights

m <- lvm(c(Y1~eta,Y2~eta,Y3~eta,eta~G+Gender))
e.lvm <- estimate(m, dW)
e.lme <- lme(value ~ time + G + Gender,
             random = ~ 1|Id,
             correlation = corSymm(),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")
e.gls <- gls(value ~ 1,#time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")

logLik(e.lvm)
logLik(e.lme)
logLik(e.gls)

test_that("mixed model UN: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                      numericDerivative = FALSE)

    expect_equal(GS.gls, res.gls)

    ## lme
    ## pb: singular information matrix
    ## GS.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
    ##                  numericDerivative = TRUE)
    ## res.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
    ##                   numericDerivative = FALSE)
})

## * latent variable model

m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
               Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)


## ** 1 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
latent(m) <- ~eta1
regression(m) <- eta1~X1+X2

e.lvm1F <- estimate(m,d)

GS.lvm1F <- dVcov2(e.lvm1F, adjust.residuals = FALSE,
                   numericDerivative = TRUE)
test_that("1 factor model: dVcov2",{    
    res.lvm1F <- dVcov2(e.lvm1F, adjust.residuals = FALSE,
                        numericDerivative = FALSE)
    expect_equal(GS.lvm1F, res.lvm1F)
    ## range(GS.lvm1F-res.lvm1F)
    ##    dfVariance(e.lvm1F)
})

## ** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e.lvm2F <- estimate(m,d)

GS.lvm2F <- dVcov2(e.lvm2F, adjust.residuals = FALSE,
                   numericDerivative = TRUE)

test_that("2 factor model: dVcov2",{
    res.lvm2F <- dVcov2(e.lvm2F, adjust.residuals = FALSE,
                        numericDerivative = FALSE)
    
    expect_equal(GS.lvm2F, res.lvm2F)
    ## range(GS.lvm2F-res.lvm2F)
})


## ** 2 factor model (covariance)

m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
## covariance(m) <- Y1 ~ Z1
latent(m) <- ~eta1+eta2

e <- estimate(m,d)

test_that("2 factor model (covariance between LV): dVcov2",{
    GS.lvm <- dVcov2(e, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    
    expect_equal(GS.lvm, res.lvm)
    ## range(GS.lvm2F-res.lvm2F)
})

## ** 2 factor model (correlation)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2+X3

e <- estimate(m,d)

test_that("2 factor model (correlation between LV): dVcov2",{
    GS.lvm <- dVcov2(e, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)
    ## range(GS.lvm2F-res.lvm2F)
})

## * constrains

m.sim <- lvm(c(Y1[mu:sigma]~X1+X2+X3,
               Y2[mu:sigma]~X1+X4+X5))
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

e.lvmC <- estimate(m.sim,d)

test_that("1 factor model: dVcov2",{    
    GS.lvmC <- dVcov2(e.lvmC, adjust.residuals = FALSE,
                       numericDerivative = TRUE)
    res.lvmC <- dVcov2(e.lvmC, adjust.residuals = FALSE,
                        numericDerivative = FALSE)
    expect_equal(GS.lvmC, res.lvmC)
})

##----------------------------------------------------------------------
### test-dVcov2.R ends here

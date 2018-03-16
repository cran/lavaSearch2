### test-Utils-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 16 2017 (10:36) 
## Version: 
## Last-Updated: mar 13 2018 (13:25) 
##           By: Brice Ozenne
##     Update #: 44
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

library(nlme)
lava.options(symbols = c("~","~~"))

context("Utils-nlme")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,Y4~1*eta,eta~G+Gender))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- lava::sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id","Gender"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]
dL$time.num <- as.numeric(dL$time)

## * Heteroschedasticity
e.gls <- nlme::gls(value ~ time + G + Gender,
                   weights = varIdent(form =~ 1|time),
                   data = dL, method = "ML")
test_that("Heteroschedasticity", {
    vec.sigma <- c(1,coef(e.gls$modelStruct$varStruct, unconstrained = FALSE))
    expect_equal(diag(vec.sigma^2 * sigma(e.gls)^2),
                 unname(getVarCov2(e.gls, cluster = "Id")$Omega))
})

## * Compound symmetry
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   data = dL,
                   method = "ML")
e.lme.bis <- nlme::lme(value ~ time + G + Gender,
                       random = ~ 1|Id,
                       correlation = corCompSymm(),
                       data = dL,
                       method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corCompSymm(form=~ 1|Id),
                   data = dL, method = "ML")

test_that("Compound symmetry", {
    expect_equal(unclass(getVarCov(e.gls)),
                 unname(getVarCov2(e.gls)$Omega))

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(getVarCov2(e.lme)$Omega))

    expect_equal(unname(getVarCov(e.lme.bis, type = "marginal", individuals = 1)[[1]]),
                 unname(getVarCov2(e.lme.bis)$Omega))
})

## * Unstructured 
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   correlation = corSymm(form =~ time.num|Id),
                   data = dL,
                   method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corSymm(form=~ time.num|Id),
                   data = dL, method = "ML")


test_that("Unstructured ", {
    expect_equal(unclass(getVarCov(e.gls)),
                 unname(getVarCov2(e.gls)$Omega))

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(getVarCov2(e.lme)$Omega))
})

## * Unstructured with weights
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   correlation = corSymm(form =~ time.num|Id),
                   weight = varIdent(form = ~ 1|time),
                   data = dL,
                   method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corSymm(form =~ time.num|Id),
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")

test_that("Unstructured with weights", {
    expect_equal(unclass(getVarCov(e.gls)),
                 unname(getVarCov2(e.gls)$Omega))

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(getVarCov2(e.lme)$Omega))
})

## * 2 random effect model (error)
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random=~1|Id/Gender,
                   data = dL,
                   method = "ML")

expect_error(getVarCov2(e.lme))

##----------------------------------------------------------------------
### test-Utils-nlme.R ends here

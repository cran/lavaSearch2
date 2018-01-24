### test-Utils-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 16 2017 (10:36) 
## Version: 
## Last-Updated: jan 19 2018 (14:45) 
##           By: Brice Ozenne
##     Update #: 34
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

context("Utils-nlme")
n <- 5e1

## * data
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,Y4~1*eta,eta~G+Gender))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- lava::sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id","Gender"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]

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

vecCoef.lme <- .coef2(e.lme)
vecCoef.lme.bis <- .coef2(e.lme.bis)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.lme.bis <- .getGroups2(e.lme.bis)
groups.gls <- .getGroups2(e.gls)

test_that("Compound symmetry", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))

    lsVcov.lme.bis <- .getVarCov2(e.lme.bis,
                                  param = vecCoef.lme.bis,
                                  attr.param = attributes(vecCoef.lme.bis),
                                  endogenous = groups.lme.bis$endogenous,
                                  name.endogenous = groups.lme.bis$name.endogenous,
                                  n.endogenous = groups.lme.bis$n.endogenous,
                                  cluster = groups.lme.bis$cluster,
                                  n.cluster = groups.lme.bis$n.cluster)

    expect_equal(unname(getVarCov(e.lme.bis, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme.bis$Omega))
})

## * Unstructured 
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   correlation = corSymm(),
                   data = dL,
                   method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corSymm(form=~ 1|Id),
                   data = dL, method = "ML")

vecCoef.lme <- .coef2(e.lme)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.gls <- .getGroups2(e.gls)

test_that("Unstructured ", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))
})

## * Unstructured with weight
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form = ~ 1|time),
                   data = dL,
                   method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corSymm(form=~ 1|Id),
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")

vecCoef.lme <- .coef2(e.lme)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.gls <- .getGroups2(e.gls)

test_that("Unstructured ", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))
})

## * 2 random effect model (error)
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random=~1|Id/Gender,
                   data = dL,
                   method = "ML")

expect_error(.getGroups2(e.lme))

## e.lme <- nlme::lme(value ~ time + G + Gender,
##              random=~1|Id,
##              correlation=corCompSymm(form = ~1|Gender),
##              data = dL,
##              method = "ML")
## incompatible

##----------------------------------------------------------------------
### test-Utils-nlme.R ends here

### test-Utils-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 16 2017 (10:36) 
## Version: 
## Last-Updated: maj  3 2018 (13:21) 
##           By: Brice Ozenne
##     Update #: 61
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

## * t.test
test_that("invariant to the order in the dataset", {
    e1.gls <- gls(Y1 ~ Gender, data = dW[order(dW$Id),],
                  weights = varIdent(form = ~1|Gender),
                  method = "ML")

    out1 <- getVarCov2(e1.gls, cluster = dW$Id)
    index.cluster <- as.numeric(names(out1$index.Omega))
    expect_true(all(diff(index.cluster)>0))

    e2.gls <- gls(Y1 ~ Gender, data = dW[order(dW$Gender),],
                  weights = varIdent(form = ~1|Gender),
                  method = "ML")
    out2 <- getVarCov2(e2.gls, cluster = dW$Id)
    index.cluster <- as.numeric(names(out2$index.Omega))
    expect_true(all(diff(index.cluster)>0))
})

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

## * PET dataset

df.PET <- data.frame("ID" = c( 925, 2020, 2059, 2051, 2072, 2156, 2159, 2072, 2020, 2051, 2231,
                              2738, 2231, 2777,  939,  539, 2738, 2777,  925, 2156, 2159, 2059), 
                     "session" = c("V", "V", "V", "V", "V", "V", "V", "C", "C", "C", "C",
                                   "C", "V", "C", "C", "V", "V", "V", "C", "C", "C", "C"), 
                     "PET" = c(-2.53, -6.74, -8.17, -2.44, -3.54, -1.27, -0.55, -0.73, -1.42,  3.35,
                               -2.11,  2.60, -4.52,  0.99, -1.02, -1.78, -5.86,  1.20, NA, NA, NA, NA)
                     )
df.PET$session.index <- as.numeric(as.factor(df.PET$session))


e.lme <- lme(PET ~ session,
             random = ~ 1 | ID,
             weights = varIdent(form=~session.index|session),
             na.action = "na.omit",
             data = df.PET)
test_that("getVarCov2 - NA", {
    expect_equal(matrix(c( 7.893839, 1.583932, 1.583932, 4.436933), 2, 2),
                 unname(getVarCov2(e.lme)$Omega), tol = 1e-6, scale = 1)
})


##----------------------------------------------------------------------
### test-Utils-nlme.R ends here

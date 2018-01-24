### test-score2-nlme.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (11:40) 
## Version: 
## last-updated: jan 19 2018 (14:45) 
##           By: Brice Ozenne
##     Update #: 103
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

library(nlme)
lava.options(symbols = c("~","~~"))

context("iid2-nlme")

## * Simulation
n <- 5e1
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- lava::sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]
dL$Z1 <- rnorm(NROW(dL))

## * Linear model
## ** Homoschedasticity (gls)
m <- lvm(c(Y1~G))
e.lvm <- estimate(m, dW)

e.gls <- nlme::gls(Y1 ~ G, data = dW, method = "ML")
allCoef.gls <- c(coef(e.gls),sigma2 = sigma(e.gls)^2)

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2.gls equivalent to score.lvm", {
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(unname(score.gls),unname(score.lvm))    
})

## ** Heteroscedasticity (gls)
m <- lvm(c(Y1~G,Y2~G,Y3~G))
e.lvm <- estimate(m, dW)

e.gls <- nlme::gls(value ~ 0+time + time:G,
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")
allCoef.gls <- .coef2(e.gls)
attributes(allCoef.gls) <- attributes(allCoef.gls)["names"]

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2.gls equivalent to score.lvm", {
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.gls.p <- score2(e.gls, p = allCoef.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    
    expect_equal(unname(score.gls[,1:6]),unname(score.lvm[,1:6]))
    expect_true(all(abs(colSums(score.gls))<1e-7))
    expect_equal(score.gls,score.gls.p)    
})

test_that("score2.gls equivalent to score.lvm (adjust residuals)", {
        score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = TRUE, indiv = TRUE)
        score.gls.p <- score2(e.gls, p = allCoef.gls, cluster = "Id", adjust.residuals = TRUE, indiv = TRUE)
        score.lvm <- score2(e.lvm, adjust.residuals = TRUE,, indiv = TRUE)

        expect_equal(unname(score.gls[,1:6]),unname(score.lvm[,1:6]))
        expect_equal(score.gls,score.gls.p)
})

## * Mixed model
## ** Compound symmetry (lme/gls)
m <- lvm(c(Y1[mu1:sigma2]~1*eta,
           Y2[mu2:sigma2]~1*eta,
           Y3[mu3:sigma2]~1*eta,
           eta~G))
e.lvm <- estimate(m, dW)

e.lme <- nlme::lme(value ~ time + G,
                   random =~1| Id,
                   data = dL, method = "ML")

e.gls <- nlme::gls(value ~ time + G,
                   correlation = corCompSymm(form = ~1| Id),
                   data = dL, method = "ML")

allCoef.lme <- .coef2(e.lme)
# .coef2(e.gls)

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lvm)))
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(allCoef.lme),as.double(coef(e.lvm)), tol = 1e-7) ## needed for CRAN
})

test_that("score2.lme/gls equivalent to score2.lvm - no adjustment", {    
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(score.lme,score.lme.p)

    score.gls <- score2(e.gls, adjust.residuals = FALSE, indiv = TRUE)
    expect_true(all(abs(colSums(score.gls))<1e-5))
    expect_equal(score.lme[,attr(allCoef.lme,"mean.coef")],
                 score.gls[,attr(allCoef.lme,"mean.coef")], tol = 1e-7)

    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(unname(score.lme),unname(score.lvm), tol = 1e-6) ## needed for CRAN
})

test_that("score2.lme/gls equivalent to score2.lvm - adjustment", {
        score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE)
        score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = TRUE, indiv = TRUE)
        expect_equal(score.lme,score.lme.p)
    
        score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE)
        expect_equal(score.lme[,attr(allCoef.lme,"mean.coef")],
                     score.gls[,attr(allCoef.lme,"mean.coef")], tol = 1e-7)

        score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE)
        expect_equal(unname(score.lme),unname(score.lvm), tol = 1e-7) ## needed for CRAN    
})

## ** Compound symmetry with weight (lme/gls)
m <- lvm(c(Y1~1*eta,
           Y2~1*eta,
           Y3~1*eta,
           eta~G))
e.lvm <- estimate(m, dW)

e.lme <- nlme::lme(value ~ time + G,
                   random =~1| Id,
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")

e.gls <- nlme::gls(value ~ time + G,
                   correlation = corCompSymm(form = ~1| Id),
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")

allCoef.lme <- .coef2(e.lme)
# .coef2(e.gls)

test_that("lme/gls equivalent to lvm", {
    ## gls not really matching
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lvm)), tol = 1e-2)

    ## lme matchin
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(allCoef.lme[c(1:5)]),as.double(coef(e.lvm)[1:5]), tol = 1e-5)
    expect_equal(as.double(allCoef.lme[8]),as.double(coef(e.lvm)[6]), tol = 1e-5)
    expect_equal(as.double(allCoef.lme[6:7]*allCoef.lme[5]),
                 as.double(coef(e.lvm)[7:8]), tol = 1e-5)
})

test_that("score2.lme/gls equivalent to score2.lvm - no adjustment", {    
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(score.lme,score.lme.p)

    score.gls <- score2(e.gls, adjust.residuals = FALSE, indiv = TRUE)
    expect_true(all(abs(colSums(score.gls))<1e-4))
    ## boxplot(score.gls) ; abline(h=0)
    ## 

    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(unname(score.lme)[,1:4],unname(score.lvm[,1:4]), tol = 1e-5)
    expect_equal(unname(score.lme)[,8],unname(score.lvm[,6]), tol = 1e-5)

    
})

test_that("score2.lme/gls equivalent to score2.lvm - adjustment", {
        score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE)
        score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = TRUE, indiv = TRUE)
        expect_equal(score.lme,score.lme.p)

        score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE)

        score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE)
        expect_equal(unname(score.lme)[,1:4],unname(score.lvm[,1:4]), tol = 1e-5)
        expect_equal(unname(score.lme)[,8],unname(score.lvm[,6]), tol = 1e-5)    
})

## ** Unstructured (lme/gls)
m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~1))
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, dW)

e.lme <- nlme::lme(value ~ time,
                         random =~1| Id,
                         correlation = corSymm(),
                         weight = varIdent(form = ~ 1|time),
                         data = dL, method = "ML")

e.gls <- nlme::gls(value ~ time,
                   correlation = corSymm(form = ~1| Id),
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")

allCoef.lme <- .coef2(e.lme)
# .coef2(e.gls)

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2.lme/gls equivalent to score2.lvm - no adjustment", {    
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = FALSE, indiv = TRUE)
    expect_equal(score.lme,score.lme.p)
    expect_true(all(abs(colSums(score.lme))<1e-3))

    score.gls <- score2(e.gls, adjust.residuals = FALSE, indiv = TRUE)
    expect_true(all(abs(colSums(score.gls))<1e-2))

    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.lme[,1:3]),unname(score.gls[,1:3]), tol = 1e-4)
    expect_equal(unname(score.lvm[,1:3]),unname(score.lme[,1:3]), tol = 1e-4)
})
#----------------------------------------------------------------------
### test-score2-nlme.R ends here

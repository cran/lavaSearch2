### test-iid2-nlme.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (12:57) 
## Version: 
## last-updated: jan 19 2018 (14:44) 
##           By: Brice Ozenne
##     Update #: 102
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

library(clubSandwich)
library(nlme)
lava.options(symbols = c("~","~~"))

context("iid2-nlme")

n <- 5e1


## * data
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,Y4~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- lava::sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]


## * Heteroscedasticity (gls)
m <- lvm(c(Y1~G,Y2~G,Y3~G,Y4~G))
e.lvm <- estimate(m, dW)

e.gls <- nlme::gls(value ~ 0+time + time:G,
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")


factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))
index.coef <- 1:length(coef(e.gls))

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("gls: HC0/HC1", {
    iid2HC0.gls <- iid2(e.gls, adjust.residuals = FALSE, data = dL, cluster = dL$Id)
    iid2HC0.lvm <- iid2(e.lvm, adjust.residuals = FALSE, cluster = dL$Id)
    expect_equal(unname(iid2HC0.gls[,index.coef]),unname(iid2HC0.lvm[,index.coef]))
    
    VsandwichHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]    
    expect_equal(unname(VsandwichHC0.gls), unname(VsandwichHC0.lvm), tolerance = 1e-10)    
    GS <- clubSandwich::vcovCR(e.gls, type = "CR0", cluster = dL$Id) * factor^2
    expect_equal(as.double(GS),as.double(VsandwichHC0.gls), tolerance = 1e-10)
    
    GS <- clubSandwich::vcovCR(e.gls, type = "CR1", cluster = dL$Id) * factor^2
    VsandwichHC1.gls <- VsandwichHC0.gls*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.gls), tolerance = 1e-10)
})

test_that("gls: HC2", {
    iid2HC2.gls <- iid2(e.gls, adjust.residuals = TRUE, data = dL, cluster = dL$Id, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, adjust.residuals = TRUE, cluster = dL$I, as.clubSandwich = 2)
    expect_equal(unname(iid2HC2.gls[,index.coef]),unname(iid2HC2.lvm[,index.coef]))
  
    VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]    
    expect_equal(unname(VsandwichHC2.gls), unname(VsandwichHC2.lvm), tolerance = 1e-10)    
    GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dL$Id) * factor^2    
    expect_equal(as.double(GS),as.double(VsandwichHC2.gls), tolerance = 1e-10) 
})

## * gls/lme/lvm - Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,Y2[mu2:sigma]~1*eta,Y3[mu3:sigma]~1*eta,Y4[mu4:sigma]~1*eta,eta~G))
e.lvm <- estimate(m, dW)

e.lme <- nlme::lme(value ~ time + G,
                   random =~1| Id,
                   data = dL, method = "ML")

e.gls <- nlme::gls(value ~ time + G,
                   correlation = corCompSymm(form = ~1| Id),
                   data = dL, method = "ML")
index.coef <- 1:length(coef(e.gls))

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})

test_that("lme: HC0/HC1", {
    iid2HC0.lme <- iid2(e.lme, data = dL, adjust.residuals = FALSE)
    iid2HC0.gls <- iid2(e.gls, data = dL, adjust.residuals = FALSE)
    iid2HC0.lvm <- iid2(e.lvm, adjust.residuals = FALSE)
    
    expect_equal(unname(iid2HC0.lme),unname(iid2HC0.lvm))
    expect_equal(unname(iid2HC0.gls[,index.coef]),unname(iid2HC0.lvm[,index.coef]))

    VsandwichHC0.lme <- crossprod(iid2HC0.lme)[index.coef,index.coef]
    VsandwichHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC0.lme), as.double(VsandwichHC0.gls))
    expect_equal(as.double(VsandwichHC0.lvm), as.double(VsandwichHC0.gls))
    GS <- clubSandwich::vcovCR(e.lme, type = "CR0", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC0.lme))
    
    GS <- clubSandwich::vcovCR(e.lme, type = "CR1", cluster = dL$Id)
    VsandwichHC1.lme <-  VsandwichHC0.lme*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.lme))    
})

test_that("lme: HC2", {
    iid2HC2.lme <- iid2(e.lme, data = dL, adjust.residuals = TRUE, as.clubSandwich = 2)
    iid2HC2.gls <- iid2(e.gls, data = dL, adjust.residuals = TRUE, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, adjust.residuals = TRUE, as.clubSandwich = 2)    
    expect_equal(unname(iid2HC2.lme),unname(iid2HC2.lvm))
    expect_equal(unname(iid2HC2.gls[,index.coef]),unname(iid2HC2.lvm[,index.coef]))

    VsandwichHC2.lme <- crossprod(iid2HC2.lme)[index.coef,index.coef]
    VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC2.lme), as.double(VsandwichHC2.gls))
    expect_equal(as.double(VsandwichHC2.lvm), as.double(VsandwichHC2.gls))
    GS <- clubSandwich::vcovCR(e.lme, type = "CR2", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC2.lme))
})

## * lme/lvm - CS with weights
m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,Y4~1*eta,eta~G))
e.lvm <- estimate(m, dW)

e.lme <- nlme::lme(value ~ time + G,
                   random =~1| Id,
                   weight = varIdent(form = ~ 1|time),
                   data = dL, method = "ML")
index.coef <- 1:length(fixef(e.lme))

## gls does not give the same likelihood
test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
#    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("lme/gls/lvm: HC0/HC1", {
    iid2HC0.lme <- iid2(e.lme, data = dL, adjust.residuals = FALSE)
    iid2HC0.lvm <- iid2(e.lvm, adjust.residuals = FALSE)    
    expect_equal(unname(iid2HC0.lme[,index.coef]),unname(iid2HC0.lvm[,index.coef]), tol = 1e-6)

    VsandwichHC0.lme <- crossprod(iid2HC0.lme)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC0.lme), as.double(VsandwichHC0.lvm), tol = 1e-7)
    GS <- clubSandwich::vcovCR(e.lme, type = "CR0", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC0.lme))
    
    GS <- clubSandwich::vcovCR(e.lme, type = "CR1", cluster = dL$Id)
    VsandwichHC1.lme <-  VsandwichHC0.lme*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.lme))    
})

test_that("lme/gls/lvm: HC2", {
    iid2HC2.lme <- iid2(e.lme, data = dL, adjust.residuals = TRUE, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, adjust.residuals = TRUE, as.clubSandwich = 2)    
    expect_equal(unname(iid2HC2.lme[,index.coef]),unname(iid2HC2.lvm[,index.coef]), tol = 1e-6)

    VsandwichHC2.lme <- crossprod(iid2HC2.lme)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC2.lme), as.double(VsandwichHC2.lvm), tol = 1e-7)
    GS <- clubSandwich::vcovCR(e.lme, type = "CR2", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC2.lme))
})

## * gls - Unstructured
m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~1))
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, dW)

dataRed <- dL[dL$time!="Y4",]
dataRed$time <- droplevels(dataRed$time)
e.gls <- nlme::gls(value ~ time,
                   correlation = corSymm(form =~ 1| Id),
                   weight = varIdent(form =~ 1|time),
                   data = dataRed, method = "ML")

e.lme <- nlme::lme(value ~ time,
                   random =~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form =~ 1|time),
                   data = dataRed, method = "ML")

index.coef <- 1:length(coef(e.gls))

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})


test_that("gls/lvm: HC0/HC1", {
    iid2HC0.gls <- iid2(e.gls, data = dataRed, adjust.residuals = FALSE)
    ## iid2HC0.lme <- iid2(e.lme, adjust.residuals = FALSE) ## not invertible
    iid2HC0.lvm <- iid2(e.lvm, adjust.residuals = FALSE)
    

    VsandwichHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]    
    expect_equal(as.double(VsandwichHC0.gls), as.double(VsandwichHC0.lvm), tol = 1e-7)
    ## scoreX <- score2(e.lvm, adjust.residuals = FALSE)
    ## colSums(scoreX)
    GS <- clubSandwich::vcovCR(e.gls, type = "CR0", cluster = dataRed$Id)
    ## VsandwichHC0.gls-GS
    ## VsandwichHC0.lvm-GS

    ## VsandwichHC0.gls-vcov(e.gls)
    ## GS-vcov(e.gls)

})

test_that("lme: HC2", {
    iid2HC2.gls <- iid2(e.gls, data = dataRed, adjust.residuals = TRUE, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, adjust.residuals = TRUE, as.clubSandwich = 2)    

    VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]    
    expect_equal(as.double(VsandwichHC2.gls), as.double(VsandwichHC2.lvm), tol = 1e-7)

    ## VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    ## GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dL$Id)
    ## expect_equal(as.double(GS),as.double(VsandwichHC2.gls))
})


#----------------------------------------------------------------------
### test-iid2-nlme.R ends here

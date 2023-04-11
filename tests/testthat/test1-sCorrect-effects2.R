### test1-sCorrect-effects2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 19 2022 (11:40) 
## Version: 
## Last-Updated: jan 19 2022 (11:57) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
## rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
context("sCorrect (effects2)")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- lava::sim(mSim, n = n, latent = FALSE)

## * effects2
m <- lvm(Y1~eta1,Y2~eta1+X2,Y3~eta1,
         Z1~eta2+X2,Z2~eta2,Z3~eta2,
         eta2~eta1+X2)

test_that("LM - effects2 correction", {

    e <- estimate(lvm(Y1~X2+Z1,Z1~X2), d)
    e2 <- estimate2(e)

    test <- effects2(e2, linfct = c("Y1~X2|total","Y1~X2|direct","Y1~X2|indirect"))
    expect_equal(summary(test, test = adjusted("none"))$table2$df,
                 c(51.82375, 47.00000, 48.04664), tol = 1e-3)
})

test_that("LVM - effects2 no correction", {

    e <- estimate(m, d)
    e0 <- estimate2(e, ssc = FALSE, df = FALSE) ## no correction

    GS <- effects(e, from = "X2", to = "Z1")
    test1 <- effects(e0, linfct = "Z1~X2")
    test3 <- effects(e0, linfct = c("Z1~X2|total","Z1~X2|direct","Z1~X2|indirect"))

    expect_equal(as.double(GS$coef[c("total","direct","indirect")]),
                 summary(test3, test = adjusted("none"))$table2[,"estimate"],
                 tol = 1e-6)
    expect_equal(as.double(summary(test1, test = adjusted("none"))$table2[1,]),
                 as.double(summary(test3, test = adjusted("none"))$table2[1,]),
                 tol = 1e-6)
    expect_equal(as.double(sqrt(diag(GS$vcov[1:3,1:3]))),
                 summary(test3, test = adjusted("none"))$table2[,"se"],
                 tol = 1e-6)


    e2 <- estimate2(e) ## correction
    test12 <- effects(e2, linfct = "Z1~X2")
    test32 <- effects(e2, linfct = c("Z1~X2|total","Z1~X2|direct","Z1~X2|indirect"))
})


##----------------------------------------------------------------------
### test1-sCorrect-effects2.R ends here

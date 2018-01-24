### test-vcov-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 13 2017 (15:52) 
## Version: 
## Last-Updated: jan 19 2018 (14:44) 
##           By: Brice Ozenne
##     Update #: 16
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

context("score2")

## * Corrected vcov
n <- 2e1

## ** linear model
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.gls <- nlme::gls(Y~X1+X2+X3, data = d, method = "ML")
e.lm <- lm(Y~X1+X2+X3, data = d)

Sigma.gls <- attr(residuals2(e.gls, data = d, cluster = 1:n, 
                             return.vcov.param = TRUE), "vcov.param")

Sigma.lm <- attr(residuals2(e.lm, data = d, 
                            return.vcov.param = TRUE), "vcov.param")
## Note: add data = d in the call even though it should not be necessary
##       however in some cases the function cannot identify the dataset in the correct environment

test_that("Corrected vcov - linear model",{
    expect_equal(unname(Sigma.gls), unname(Sigma.lm))
})


##----------------------------------------------------------------------
### test-vcov-nlme.R ends here

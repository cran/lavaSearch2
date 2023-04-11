### test1-sCorrect-validObjects.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  6 2018 (10:42) 
## Version: 
## Last-Updated: Jan 11 2022 (09:07) 
##           By: Brice Ozenne
##     Update #: 68
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
library(data.table)
lava.options(symbols = c("~","~~"))
context("sCorrect (warnings and errors for invalid objects/arguments)")

## * Simulation
n <- 100
m.sim <- lvm(Y~X1+X2,G~1)
categorical(m.sim,K=3,label=c("a","b","c")) <- ~G+X2
set.seed(10)
d <- lava::sim(m.sim,n,latent=FALSE)

## * estimate2 for lvm objects

## ** error for multigroup lvm
## check in sCorrect.R 
suppressWarnings(e <- estimate(list(lvm(Y~X1),lvm(Y~X1),lvm(Y~X1)), data = split(d,d$G)))
test_that("error for multigroup models", {
    expect_error(estimate2(e))
})

## ** error for tobit lvm
## check in sCorrect.R
e <- estimate(lvm(G~X1), data = d)
test_that("error for tobit models", {
    expect_error(estimate2(e))
})

## ** error for lvm with transform variables
## check in sCorrect.R
m <- lvm(Y~X1)
transform(m,Id~X1) <- function(x){1:NROW(x)}
d.tempo <- lava::sim(m, n)
e <- estimate(m, data = d.tempo)
test_that("error when using transform", {
    expect_error(estimate2(e))
})

## * estimate2 with data.table
dt <- as.data.table(d)
e <- estimate(lvm(Y~X1+X2+G), data = dt)
test_that("ok for data.table objects", {
    expect_true(inherits(estimate2(e), "lvmfit2"))
})

##----------------------------------------------------------------------
### test1-sCorrect-validObjects.R ends here


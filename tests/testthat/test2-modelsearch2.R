### test2-modelsearch2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2018 (11:45) 
## Version: 
## Last-Updated: jan 22 2018 (12:11) 
##           By: Brice Ozenne
##     Update #: 6
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

lava.options(symbols = c("~","~~"))

context("lavaSearch2")

## * example with one additional link
n <- 100
m.sim <- lvm(Y~E+0*X1)
m <- lvm(Y~E)
addvar(m) <- ~X1

set.seed(10)
df.sim <- sim(m.sim, n=100, latent = FALSE)
e.base <- estimate(m, data = df.sim)

test_that("Score 1 link",{
    GS.score <- modelsearch(e.base, silent = TRUE)
    index.coef <- which(GS.score$res[,"Index"]=="Y~X1")

    test.score <- summary(modelsearch2(e.base, statistic = "score", method.p.adjust = "holm", trace = 0), display = FALSE)

    
    expect_equivalent(GS.score$test[index.coef,"Test Statistic"],test.score$data[1,"statistic"])
    expect_equivalent(GS.score$test[index.coef,"P-value"],test.score$data[1,"p.value"])
})

test_that("Wald 1 link",{
    ## no adjustement
    e.GS <-  estimate(lvm(Y~E+X1), data = df.sim)
    test.GS <- lTest(e.GS, adjust.residuals = FALSE)
    
    test.Wald <- summary(modelsearch2(e.base, link = "Y~X1", df = FALSE, adjust.residuals = FALSE,
                                       statistic = "Wald", method.p.adjust = "max", trace = 0), display = FALSE)

    expect_equivalent(abs(test.GS["Y~X1","statistic"]), test.Wald$data[1,"statistic"])
    expect_equal(test.GS["Y~X1","p-value"], test.Wald$data[1,"p.value"], tolerance = 1e-4)

    ## with adjustment
    e.GS <-  estimate(lvm(Y~E+X1), data = df.sim)
    test.GS <- lTest(e.GS, adjust.residuals = TRUE)
    
    test.Wald <- summary(modelsearch2(e.base, link = "Y~X1", df = TRUE, adjust.residuals = TRUE,
                                       statistic = "Wald", method.p.adjust = "max", trace = 0), display = FALSE)

    expect_equivalent(abs(test.GS["Y~X1","statistic"]), test.Wald$data[1,"statistic"])
    expect_equal(test.GS["Y~X1","p-value"], test.Wald$data[1,"p.value"], tolerance = 1e-4)
})

##----------------------------------------------------------------------
### test2-modelsearch2.R ends here

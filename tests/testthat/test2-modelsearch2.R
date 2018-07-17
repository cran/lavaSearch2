### test2-modelsearch2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2018 (11:45) 
## Version: 
## Last-Updated: maj  2 2018 (09:55) 
##           By: Brice Ozenne
##     Update #: 18
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

lava.options(symbols = c("~","~~"))

context("lavaSearch2")

## * example with one additional link
n <- 100
m.sim <- lvm(Y~E+0*X1)
m <- lvm(Y~E)
addvar(m) <- ~X1

set.seed(10)
df.sim <- lava::sim(m.sim, n=100, latent = FALSE)
e.base <- estimate(m, data = df.sim)

test_that("Score 1 link",{
    GS.score <- modelsearch(e.base, silent = TRUE)
    index.coef <- which(GS.score$res[,"Index"]=="Y~X1")

    test.score <- summary(modelsearch2(e.base, statistic = "score", method.p.adjust = "holm", trace = 0), print = FALSE)

    
    expect_equivalent(GS.score$test[index.coef,"Test Statistic"],test.score$data[1,"statistic"])
    expect_equivalent(GS.score$test[index.coef,"P-value"],test.score$data[1,"p.value"])
})

test_that("Wald 1 link",{
    ## no adjustement
    e.GS <-  estimate(lvm(Y~E+X1), data = df.sim)
    test.GS <- compare2(e.GS, par = "Y~X1" , bias.correct = FALSE, as.lava = FALSE)
    
    test.Wald <- summary(modelsearch2(e.base, link = "Y~X1", df = FALSE, bias.correct = FALSE,
                                       statistic = "Wald", method.p.adjust = "max", trace = 0), print = FALSE)

    expect_equivalent(abs(test.GS[1,"statistic"]), test.Wald$data[1,"statistic"])
    expect_equal(test.GS[1,"p-value"], test.Wald$data[1,"p.value"], tolerance = 1e-5)

    ## with adjustment
    e.GS <-  estimate(lvm(Y~E+X1), data = df.sim)
    test.GS <- compare2(e.GS, par  = "Y~X1", bias.correct = TRUE, as.lava = FALSE)
    
    test.Wald <- summary(modelsearch2(e.base, link = "Y~X1", df = TRUE, bias.correct = TRUE,
                                       statistic = "Wald", method.p.adjust = "max", trace = 0), print = FALSE)

    expect_equivalent(abs(test.GS[1,"statistic"]), test.Wald$data[1,"statistic"])
    expect_equal(test.GS[1,"p-value"], test.Wald$data[1,"p.value"], tolerance = 1e-5)
})

##----------------------------------------------------------------------
### test2-modelsearch2.R ends here

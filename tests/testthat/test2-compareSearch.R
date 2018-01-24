### test-compareSearch.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (09:25) 
## Version: 
## last-updated: jan 18 2018 (17:51) 
##           By: Brice Ozenne
##     Update #: 19
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * header
if(TRUE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

library(mvtnorm)
lava.options(symbols = c("~","~~"))

context("compareSearch")

## * example
n <- 100
m.sim <- lvm(Y~E+X1+0.5*X2+0.1*X3,~Z1+Z2+Z3)
m.base <- lvm(Y~E,~X1+X2+X3+Z1+Z2+Z3)

set.seed(10)
df.sim <- sim(m.sim, n=100, latent = FALSE)
e.base <- estimate(m.base, data = df.sim)

resCompare <- compareSearch(e.base, 
                            method.p.adjust = c("none","bonferroni","fdr","max"),
                            statistic = c("score","Wald"),trace = 0)

resCompare <- compareSearch(e.base, link = c("Y~X1","Y~X3","Y~Z1"),
                            method.p.adjust = c("none","bonferroni","fdr"),
                            statistic = c("score"),trace = 0)

## * example with error (wrong link)
resCompare <- compareSearch(e.base, link = c("Y~G","Y~A"),
                            method.p.adjust = c("none","bonferroni","fdr","max"),
                            statistic = c("score","Wald"),trace = 0)
##----------------------------------------------------------------------
### test-compareSearch.R ends here

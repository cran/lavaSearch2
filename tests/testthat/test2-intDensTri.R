### test-intDensTri.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 31 2017 (16:32) 
## Version: 
## last-updated: feb  6 2018 (09:29) 
##           By: Brice Ozenne
##     Update #: 14
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
    rm(list = ls(all.names = TRUE))
    library(testthat)
    library(lavaSearch2)
}

library(mvtnorm)
lava.options(symbols = c("~","~~"))

context("intDensTri")

## * tests
# {{{ around 0
test_that("Integrate standard gaussian density (2D)", {
    p <- 2
    Sigma <- diag(p)
    mu <- rep(0, p)

    for(n in c(5,10,20,50,100)){ # n <- 5
        res <- intDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0)
        expect_equal(res$value,
                     1/2,
                     tol = 1e-6)
    }
})

test_that("Integrate standard gaussian density (3D)", {
    p <- 3
    Sigma <- diag(p)
    mu <- rep(0, p)

    for(n in c(5,10,20,50,100)){
        res <- intDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0, z.max=10)
        expect_equal(res$value,
                     1/2,
                     tol = 1e-6)
    }
})
# }}}

# {{{ far from 0
test_that("Integrate standard gaussian density (2D)", {
    p <- 2
    Sigma <- diag(p)
    mu <- c(10,0)

    for(n in c(5,10,20,50,100)){
        res <- intDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0)
        expect_equal(res$value,
                     1,
                     tol = 1e-6)
    }
})

test_that("Integrate standard gaussian density (3D)", {
    p <- 3
    Sigma <- diag(p)
    mu <- c(10,0,0)

    for(n in c(5,10,20,50,100)){
        res <- intDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0, z.max=10)
        expect_equal(res$value,
                     1,
                     tol = 1e-6)
    }
})
# }}}


#----------------------------------------------------------------------
### test-intDensTri.R ends here

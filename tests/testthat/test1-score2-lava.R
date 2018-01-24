### test-score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 13 2017 (11:28) 
## Version: 
## last-updated: jan 19 2018 (15:18) 
##           By: Brice Ozenne
##     Update #: 189
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

library(lava.tobit)
lava.options(symbols = c("~","~~"))

context("score2-lava")
n <- 5e1

## * linear regression
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("linear regression (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("linear regression (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("linear regression (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)
})

test_that("linear regression: constrains",{
    m <- lvm(Y[0:2]~X1+1*X2)
    e <- estimate(m, d)
    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
    
    m <- lvm(Y~beta*X1+beta*X2)
    e <- estimate(m, d)
    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
})


## * multiple linear regression
## ** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("multiple linear regression (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("multiple linear regressions: constrains",{
    m <- lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
})

## ** with covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("multiple linear regression, covariance link (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)    
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression, covariance link (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression, covariance link (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)
})

## * latent variable model
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## ** factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("factor model (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
     
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})
test_that("factor model (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("factor model (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("factor model: fixed coefficients",{
    m <- lvm(Y1~1*eta+1*X2,Y2~1*eta,Y3~1*eta)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
})

test_that("factor model: constrains",{
    m <- lvm(Y1~1*eta+X2,Y2~lambda*eta+X2,Y3~lambda*eta,eta ~ beta*X2+beta*X1)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
})


## ** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("2 factor model: constrains",{
    m <- lvm(Y1~1*eta1+X2,Y2~lambda*eta1+X2,Y3~lambda*eta1,eta1 ~ beta*X2+beta*X1,
             Z1~0+eta2,Z2~lambda*eta2,Z3~eta2)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
})

## ** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
## covariance(m) <- Y1 ~ Z1
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model, covariance (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariance (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariance (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

## ** 2 factor model (correlation LV)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2+X3

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- FALSE

test_that("2 factor model, correlation LV (at ML)",{
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, correlation LV (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariancecorrelation LV (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

## * leverage adjusted score for lvm models
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## ** linear regression
e0 <- estimate(lvm(Y1~X1),d)
s0 <- score2(e0, adjust.residuals = TRUE)

test_that("score2.lvm vs score2.lm (adj)",{
    sGS <- score2(lm(Y1~X1, data = d), adjust.residuals = TRUE)
    expect_equal(unname(s0),unname(sGS))
})

## ** multiple linear regression
e1 <- estimate(lvm(Y1~X1,Y2~X2+X3,Y3~1),d)

test_that("score2.lvm(adj): univariate vs. multiple univariate",{
    s1 <- score2(e1, adjust.residuals = TRUE)
    expect_equal(s1[,c("Y1","Y1~X1","Y1~~Y1")],
                 s0[,c("Y1","Y1~X1","Y1~~Y1")])
})



## ** lvm model with leverage adjusted residuals
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1

e <- estimate(m,d)

test_that("",{
    r2 <- score2(e)
})


## * error for tobit and multigroup lvm
## ** multigroup model
m.sim <- lvm(Y~X1+X2,G~1)
categorical(m.sim,K=2,label=c("a","b")) <- ~G+X2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(list(m.sim,m.sim),data = split(d,d$G))

test_that("error for multigroup models", {
    expect_error(score2(e))
})

## ** model with binary endogenous variables
m.sim <- lvm(Y~X1)
categorical(m.sim,K=2,labels = c("a","b")) <- ~Y
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(lvm(Y~X1),data = d)

test_that("error for tobit models", {
    expect_error(score2(e))
})

##----------------------------------------------------------------------
### test-score2.R ends here


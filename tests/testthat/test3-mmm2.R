### test-mmm2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (15:22) 
## Version: 
## Last-Updated: jan 18 2018 (17:50) 
##           By: Brice Ozenne
##     Update #: 54
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
if(TRUE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

library(multcomp)
library(sandwich)
lava.options(symbols = c("~","~~"))

context("multcomp - mmm")

## * simulation
mSim <- lvm(c(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10)~ beta * eta, E ~ 1)
latent(mSim) <- "eta"
set.seed(10)
n <- 1e2

df.data <- sim(mSim, n, latent = FALSE, p = c(beta = 1))

## * linear regression
name.Y <- setdiff(endogenous(mSim),"E")
n.Y <- length(name.Y)

ls.formula <- lapply(paste0(name.Y,"~","E"),as.formula)
ls.lm <- lapply(ls.formula, lm, data = df.data)
names(ls.lm) <- name.Y



test_that("mmm2 vs mmm", {
    ## rescaling using std dev
    class(ls.lm) <- "mmm"
    e.glht <- glht(ls.lm, mlf("E = 0"))

    class(ls.lm) <- "mmm2"
    e.glht2 <- glht(ls.lm, mlf2("E = 0"),
                    adjust.residuals = FALSE, robust = FALSE)

    expect_equal(e.glht$vcov, n/(n-2)*e.glht2$vcov)
    e.glht$vcov <- NULL
    e.glht2$vcov <- NULL
    ##e.glht$df <- 0    
    e.glht2$df <- 0
    e.glht2$model <- NULL
    e.glht$model <- NULL
    expect_equal(e.glht,e.glht2)

    ## no rescaling
    class(ls.lm) <- "mmm"
    e.glht <- glht(ls.lm, mlf("E = 0"),
                   vcov = sandwich)
    
    class(ls.lm) <- "mmm2"    
    e.glht2 <- glht(ls.lm, mlf2("E = 0"),
                    adjust.residuals = FALSE, robust = TRUE)
    
    e.glht2$df <- 0
    e.glht2$model <- NULL
    e.glht$model <- NULL
    dimnames(e.glht$vcov) <- NULL
    dimnames(e.glht2$vcov) <- NULL
    expect_equal(e.glht,e.glht2)

    ## same p.values
    system.time(
        res.GS <- summary(e.glht)
    )

    iid.tempo <- do.call(cbind,lapply(ls.lm, iid)) %*% t(e.glht$linfct)
    beta <- unlist(lapply(ls.lm, coef)) %*% t(e.glht$linfct)
    beta.var <- diag(crossprod(iid.tempo))
    z.value <- beta/sqrt(beta.var)
    system.time(
        res.Search <- calcDistMaxIntegral(as.vector(z.value),
                                          iid = iid.tempo, quantile.compute = FALSE,
                                          df = NULL, trace = FALSE, alpha = 0.05)
    )
    expect_equal(as.double(res.Search$p.adjust),
                 as.double(res.GS$test$pvalues),
                 tol = attr(res.GS$test$pvalues, "error")
                 )
})


## * lvm
ls.lvm <- list(Y1 = estimate(lvm(Y1~E), data = df.data),
               Y2 = estimate(lvm(Y2~E), data = df.data),
               Y3 = estimate(lvm(Y3~E), data = df.data),
               Y4 = estimate(lvm(Y4~E), data = df.data),
               Y5 = estimate(lvm(Y5~E), data = df.data),
               Y6 = estimate(lvm(Y6~E), data = df.data),
               Y7 = estimate(lvm(Y7~E), data = df.data),
               Y8 = estimate(lvm(Y8~E), data = df.data),
               Y9 = estimate(lvm(Y9~E), data = df.data),
               Y10 = estimate(lvm(Y10~E), data = df.data)
               )

ls.lm <- list(Y1 = lm(Y1~E, data = df.data),
              Y2 = lm(Y2~E, data = df.data),
              Y3 = lm(Y3~E, data = df.data),
              Y4 = lm(Y4~E, data = df.data),
              Y5 = lm(Y5~E, data = df.data),
              Y6 = lm(Y6~E, data = df.data),
              Y7 = lm(Y7~E, data = df.data),
              Y8 = lm(Y8~E, data = df.data),
              Y9 = lm(Y9~E, data = df.data),
              Y10 = lm(Y10~E, data = df.data)
              )
class(ls.lm) <- "mmm"

test_that("ls.lvmfit vs mmm", {

    ##
    C <- createContrast(ls.lvm, var.test = "E")
    lvm2.glht <- glht2(ls.lvm, linfct = C,
                       adjust.residuals = FALSE, robust = TRUE)
    lvm2.sglht <- summary(lvm2.glht)    

    ##
    class(ls.lvm) <- "ls.lvmfit"
    lvm.coef <- names(unlist(lapply(ls.lvm,coef)))
    target.coef <- grep("E",lvm.coef, value = TRUE)
    n.coef <- length(lvm.coef)
    n.target <- length(target.coef)
    lvm.C <- matrix(0, n.target, n.coef, dimnames = list(target.coef, lvm.coef) )
    diag(lvm.C[target.coef,target.coef]) <- 1

    lvm.glht <- glht(ls.lvm, linfct = lvm.C)
    lvm.glht$vcov <- vcov(ls.lvm, return.null = FALSE,
                          adjust.residuals = FALSE, robust = TRUE)
    lvm.glht$df <- NROW(df.data)
    lvm.sglht <- summary(lvm.glht)

    ## mmm
    lm.coef <- names(unlist(lapply(ls.lm,coef)))
    target.coef <- grep("E",lm.coef, value = TRUE)
    n.target <- length(target.coef)
    n.coef <- length(lm.coef)
    lm.C <- matrix(0, nrow = n.target, ncol = n.coef,
                   dimnames = list(target.coef, lm.coef) )
    diag(lm.C[target.coef,target.coef]) <- 1
    
    lm.glht <- glht(ls.lm, linfct = lm.C, vcov = sandwich)
    lm.sglht <- summary(lm.glht)

    ## compare
    expect_equal(as.numeric(lvm2.sglht$test$coefficients),
                 as.numeric(lvm.sglht$test$coefficients))

    expect_equal(as.numeric(lvm2.sglht$test$sigma),
                 as.numeric(lvm.sglht$test$sigma))

    expect_equal(as.numeric(lvm2.sglht$test$pvalues),
                 as.numeric(lvm.sglht$test$pvalues),
                 tol = attr(lm.sglht$test$pvalues,"error")*10)
    
    expect_equal(as.double(lvm.sglht$test$coefficients),
                 as.double(lm.sglht$test$coefficients))

    expect_equal(as.double(lvm.sglht$test$sigma),
                 as.double(lm.sglht$test$sigma))

    expect_equal(as.double(lvm.sglht$test$pvalues),
                 as.double(lm.sglht$test$pvalues),
                 tol = attr(lm.sglht$test$pvalues,"error")*10)
})

##----------------------------------------------------------------------
### test-mmm2.R ends here

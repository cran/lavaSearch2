### test-lTest.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: jan 19 2018 (16:02) 
##           By: Brice Ozenne
##     Update #: 153
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
library(lme4)
library(lmerTest)
library(pbkrtest)
lava.options(symbols = c("~","~~"))

context("lTest")

## * Simulation
n <- 5e1
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender,X1~1,X2~1))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,id.vars = c("G","Id","Gender","X1","X2"), variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]

## * linear regression
## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## using the t test function
e.ttest <- t.test(dW$Y1,dW$Y2)
e.ttest$parameter

## by hand
sX1 <- var(dW$Y1)/n
sX2 <- var(dW$Y2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))

df-e.ttest$parameter

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2), data = dW)
e.lm <- lm(Y1~X1+X2, data = dW)
e.gls <- nlme::gls(Y1~X1+X2, data = dW, method = "ML")

## vcov(e.lvm)

### *** clubSandwich
cS.vcov <- vcovCR(e.lm, type = "CR0", cluster = dW$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(dW))
cS.df
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p

### *** lTest
test_that("linear regression: df",{

    df.lvm <- lTest(e.lvm, adjust.residuals = FALSE, Ftest = FALSE)
    df.lm <- lTest(e.lm, adjust.residuals = FALSE, Ftest = FALSE)
    df.gls <- lTest(e.gls, cluster = 1:n, adjust.residuals = FALSE, Ftest = FALSE)
    
    ## test equivalence
    expect_equivalent(df.lvm,df.gls)
    expect_equivalent(df.lvm,df.lm)

    ## test value
    n.param <- length(coef(e.lm))
    GS <- c(rep(NROW(dW),n.param), NROW(dW)/4)
    expect_equal(df.lm$df,GS)
})

test_that("linear regression: df adjusted",{
    df.adj.lm <- lTest(e.lm, adjust.residuals = TRUE, Ftest = FALSE)
    df.adj.lvm <- lTest(e.lvm, adjust.residuals = TRUE, Ftest = FALSE)
    df.adj.gls <- lTest(e.gls, cluster = 1:n, adjust.residuals = TRUE, Ftest = FALSE)
    
    ## test equivalence
    expect_equivalent(df.adj.lvm,df.adj.gls)
    expect_equivalent(df.adj.lvm,df.adj.lm)
    
    ## test value
    n.param <- length(coef(e.lm))
    GS <- c(rep(NROW(dW)-n.param,n.param), (NROW(dW)-n.param)/4)
    expect_equal(df.adj.lm$df,GS)    
    
})

## * mixed model

## ** Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender)) 
e.lvm <- estimate(m, dW)
## lTest(e.lvm)

e.lmer <- lme4::lmer(value ~ time + G + Gender + (1|Id),
               data = dL, REML = FALSE)

e.lme <- nlme::lme(value ~ time + G + Gender, random = ~ 1|Id, data = dL, method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corCompSymm(form=~ 1|Id),
                   data = dL, method = "ML")

## *** clubSandwich - bug
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dL$Id)
## strange that same type of coef have very different degrees of freedom


## *** lava - ok
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed model: df",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Satterthwaite")$df

    df1.lvm <- lTest(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = FALSE)
    df2.lvm <- lTest(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    expect_equal(df1.lvm,df2.lvm)
    expect_equal(as.double(GS),
                 as.double(df1.lvm[1:5,"df"]), tol = 1e-4) ## needed for CRAN

    df1.lme <- lTest(e.lme, adjust.residuals = FALSE,
                          numericDerivative = FALSE)
    df2.lme <- lTest(e.lme, adjust.residuals = FALSE,
                          numericDerivative = TRUE)
    expect_equal(df1.lme, df2.lme)
    expect_equal(unname(GS), df1.lme[1:5,"df"], tol = 1e-4) ## needed for CRAN

    df1.gls <- lTest(e.gls, adjust.residuals = FALSE,
                          numericDerivative = FALSE)
    df2.gls <- lTest(e.gls, adjust.residuals = FALSE,
                          numericDerivative = TRUE)
    expect_equal(df1.gls, df2.gls) 
    expect_equal(unname(GS), df1.gls[1:5,"df"], tol = 1e-4)

})

test_that("mixed model: df adjusted",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Kenward-Roger")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Kenward-Roger")$df

    ## get_Lb_ddf(e.lmer, c(0,1,0,0,0))
    ## get_Lb_ddf(e.lmer, c(0,0,0,1,0))
    df.adj.lvm <- lTest(e.lvm, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    df.adj.lvm

    df.adj.lme <- lTest(e.lme, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    df.adj.lme

    df.adj.gls <- lTest(e.gls, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    df.adj.gls
})

## ** Unstructured with weights

m <- lvm(c(Y1~eta,Y2~eta,Y3~eta,eta~G+Gender))
e.lvm <- estimate(m, dW)
e.lme <- nlme::lme(value ~ time + G + Gender,
                   random = ~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form = ~1|time),
                   data = dL, method = "ML")
e.gls <- nlme::gls(value ~ time + G + Gender,
                   correlation = corSymm(form=~ 1|Id),
                   weight = varIdent(form = ~1|time),
                   data = dL, method = "ML")

logLik(e.lvm)
logLik(e.lme)
logLik(e.gls)

test_that("UN mixed model: df",{
    ## singular information matrix
    ## df.adj.lme <- lTest(e.lme,
    ##                          robust = FALSE, adjust.residuals = FALSE)
    skip_on_cran()
    
    system.time(
        df1.gls <- lTest(e.gls, adjust.residuals = FALSE,
                         numericDerivative = TRUE)
    )
    system.time(
        df2.gls <- lTest(e.gls, adjust.residuals = FALSE,
                         numericDerivative = FALSE)
    )    
    expect_equal(df1.gls,df2.gls)
    system.time(
        df2.adj.gls <- lTest(e.gls, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    )
    df2.adj.gls

    system.time(
        df1.lvm <- lTest(e.lvm, adjust.residuals = FALSE,
                         numericDerivative = TRUE)
    )
    system.time(
        df2.lvm <- lTest(e.lvm, adjust.residuals = FALSE,
                         numericDerivative = FALSE)
    )
    expect_equal(df1.lvm,df2.lvm)
    system.time(
        df2.adj.lvm <- lTest(e.lvm, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    )
    df2.adj.lvm
})


#----------------------------------------------------------------------
### test-lTest.R ends here


### test-compare2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: feb  5 2018 (13:48) 
##           By: Brice Ozenne
##     Update #: 191
#----------------------------------------------------------------------
## 
### Commentary: 
## Test battery:
## - linear regression:
##   compare2.lvm and compare2.lm
##   * check whether df and sigma match with what is expected
##     with and without small sample correction
##
## - mixed model (compound symmetry)
##   compare2.lvm, compare2.gls, and compare2.lme
##   * compare Satterthwaite correction with lmerTest (Wald test)
##   * compare with lava::compare (F-test)
##   * small sample correction: assess internal consistency between compare2.lvm, compare2.gls, and compare2.lme
##
## - mixed model (unstructured)
##   compare2.lvm, compare2.gls, and compare2.lme
##   * Satterthwaite correction: consistency between using or not numerical derivative
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
.coef2 <- lavaSearch2:::.coef2

context("compare2")

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

### *** compare2
test_that("linear regression: df",{
    name.param <- names(coef(e.lvm))
    n.param <- length(name.param)        
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:n.param,]

    name.param <- names(.coef2(e.gls))
    n.param <- length(name.param)        
    df.lm <- compare2(e.lm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:n.param,]

    df.gls <- compare2(e.gls, par = name.param, cluster = 1:n, bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    
    ## test equivalence
    expect_equivalent(df.lvm,df.gls)
    expect_equivalent(df.lvm,df.lm)

    ## test value
    n.param <- length(coef(e.lm))
    df.GS <- c(rep(NROW(dW),n.param), NROW(dW)/4)
    expect_equal(df.lm$df, df.GS)

    sigma2 <- coef(e.lvm)["Y1~~Y1"]
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/e.lvm$data$n))
    expect_equal(df.lm$std, unname(std.GS))
})

test_that("linear regression: df adjusted",{
    name.param <- names(coef(e.lvm))
    n.param <- length(name.param)        
    df.adj.lvm <- compare2(e.lvm, par = name.param,
                           bias.correct = TRUE, as.lava = FALSE)[1:n.param,]

    name.param <- names(.coef2(e.gls))
    n.param <- length(name.param)        
    df.adj.lm <- compare2(e.lm, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:n.param,]

    df.adj.gls <- compare2(e.gls, par = name.param, cluster = 1:n, bias.correct = TRUE, as.lava = FALSE)[1:n.param,]
    
    ## test equivalence
    expect_equivalent(df.adj.lvm,df.adj.gls)
    expect_equivalent(df.adj.lvm,df.adj.lm)
    
    ## test value
    n.param <- length(coef(e.lm))
    GS <- c(rep(NROW(dW)-n.param,n.param), (NROW(dW)-n.param)/4)
    expect_equal(df.adj.lm$df,GS)    

    sigma2 <- (e.lvm$data$n+n.param)/e.lvm$data$n * coef(e.lvm)["Y1~~Y1"]
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/(e.lvm$data$n-n.param)))
    expect_equal(df.adj.lm$std, unname(std.GS))
})

## * mixed model

## ** Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender)) 
e.lvm <- estimate(m, dW)
## compare2(e.lvm)

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

test_that("mixed model: Satterthwaite correction (Wald)",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Satterthwaite")

    name.param <- names(coef(e.lvm))
    n.param <- length(name.param)        
    df1.lvm <- compare2(e.lvm, par = name.param, numeric.derivative = FALSE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    df2.lvm <- compare2(e.lvm, par = name.param, numeric.derivative = TRUE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    expect_equal(df1.lvm,df2.lvm)
    expect_equal(as.double(GS$df),
                 as.double(df1.lvm[1:5,"df"]), tol = 1e-4) ## needed for CRAN
    expect_equal(as.double(GS$tvalue),
                 as.double(df1.lvm[1:5,"statistic"]), tol = 1e-8) ## needed for CRAN

    name.param <- names(.coef2(e.lme))
    n.param <- length(name.param)        
    df1.lme <- compare2(e.lme, par = name.param, numeric.derivative = FALSE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    df2.lme <- compare2(e.lme, par = name.param, numeric.derivative = TRUE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    expect_equal(df1.lme, df2.lme)
    expect_equivalent(df1.lme, df1.lvm)

    name.param <- names(.coef2(e.gls))
    n.param <- length(name.param)        
    df1.gls <- compare2(e.gls, par = name.param, numeric.derivative = FALSE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    df2.gls <- compare2(e.gls, par = name.param, numeric.derivative = TRUE,
                        bias.correct = FALSE, as.lava = FALSE)[1:n.param,]
    expect_equal(df1.gls, df2.gls)
    expect_equal(unlist(df1.gls[1:5,]), unlist(df1.lvm[1:5,]), tolerance = 1e-6)
})

test_that("mixed model: Satterthwaite correction (F-test)",{
    GS <- lava::compare(e.lvm, par = c("eta~G"))
    e.lava2 <- compare2(e.lvm, par = c("eta~G"), bias.correct = FALSE, as.lava = FALSE)
    expect_equal(e.lava2["global","statistic"], as.double(GS$statistic))

    GS <- lava::compare(e.lvm, par = c("Y2","Y3"))
    e.lava2 <- compare2(e.lvm, par = c("Y2","Y3"), bias.correct = FALSE, as.lava = FALSE)
    expect_equal(2*e.lava2["global","statistic"], as.double(GS$statistic))
})


test_that("mixed model: KR-like correction",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Kenward-Roger")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Kenward-Roger")

    ## get_Lb_ddf(e.lmer, c(0,1,0,0,0))
    ## get_Lb_ddf(e.lmer, c(0,0,0,1,0))
    name.param <- names(coef(e.lvm))
    n.param <- length(name.param)        
    df.adj.lvm <- compare2(e.lvm, par = name.param, numeric.derivative = FALSE,
                           bias.correct = TRUE, as.lava = FALSE)

    name.param <- names(.coef2(e.lme))
    n.param <- length(name.param)        
    df.adj.lme <- compare2(e.lme, par = name.param, numeric.derivative = FALSE,
                           bias.correct = TRUE, as.lava = FALSE)
    expect_equivalent(df.adj.lvm, df.adj.lme)

    name.param <- names(.coef2(e.gls))
    n.param <- length(name.param)        
    df.adj.gls <- compare2(e.gls, par = name.param, numeric.derivative = FALSE,
                           bias.correct = TRUE, as.lava = FALSE)
    expect_equal(unlist(df.adj.gls[1:5,]), unlist(df.adj.lvm[1:5,]), tolerance = 1e-6)
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
    ## df.adj.lme <- compare2(e.lme,
    ##                          robust = FALSE, bias.correct = FALSE)
    skip_on_cran()
    
    name.param <- names(coef(e.lvm))
    n.param <- length(name.param)        
    df1.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE,
                        numeric.derivative = FALSE, as.lava = FALSE)
    df2.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE,
                        numeric.derivative = TRUE, as.lava = FALSE)
    expect_equal(df1.lvm, df2.lvm)

    name.param <- names(.coef2(e.gls))
    n.param <- length(name.param)        
    df1.gls <- compare2(e.gls, par = name.param, bias.correct = FALSE,
                        numeric.derivative = TRUE, as.lava = FALSE)
    df2.gls <- compare2(e.gls, par = name.param, bias.correct = FALSE,
                        numeric.derivative = FALSE, as.lava = FALSE)
    expect_equal(df1.gls,df2.gls)


})


#----------------------------------------------------------------------
### test-compare2.R ends here


### test-compare2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: mar 13 2018 (16:36) 
##           By: Brice Ozenne
##     Update #: 212
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
rm(list = ls())
if(FALSE){ ## already called in test-all.R
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
.coef2.gls <- lavaSearch2:::.coef2.gls
.coef2.lme <- lavaSearch2:::.coef2.lme
context("compare2")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]

## * linear regression [lm,gls,lvm]

e.lvm <- estimate(lvm(Y1~X1+X2), data = d)
e.lm <- lm(Y1~X1+X2, data = d)
e.gls <- nlme::gls(Y1~X1+X2, data = d, method = "ML")

## vcov(e.lvm)

### ** clubSandwich
cS.vcov <- vcovCR(e.lm, type = "CR0", cluster = d$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(d))
cS.df
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p

### ** compare2
test_that("linear regression: df",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    
    name.param <- names(.coef2(e.gls))
    df.lm <- compare2(e.lm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    df.gls <- compare2(e.gls, par = name.param, cluster = 1:n, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    
    ## test equivalence
    expect_equivalent(df.lvm,df.gls)
    expect_equivalent(df.lvm,df.lm)

    ## test value
    n.param <- length(coef(e.lm))
    df.GS <- c(rep(n,n.param), n/4)
    expect_equal(df.lm$df, df.GS)

    sigma2 <- coef(e.lvm)["Y1~~Y1"]
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/e.lvm$data$n))
    expect_equal(df.lm$std, unname(std.GS))
})

test_that("linear regression: df adjusted",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]
    
    name.param <- names(.coef2(e.gls))
    df.lm <- compare2(e.lm, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]
    df.gls <- compare2(e.gls, par = name.param, cluster = 1:n, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]

    ## test equivalence
    expect_equivalent(df.lvm,df.gls)
    expect_equivalent(df.lvm,df.lm)

    ## test value
    n.param <- length(coef(e.lm))
    df.GS <- c(rep(n-n.param,n.param), (n-n.param)/4)
    expect_equal(df.lm$df, df.GS)

    sigma2 <- sigma(e.lm)^2
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/(n-n.param)))
    expect_equal(df.lm$std, unname(std.GS), tolerance = 1e-7)
})

## * multiple linear regression [lvm,gls]
## ** model fit
ls.lm <- list(lm(Y1~X1,d),lm(Y2~X2,d),lm(Y3~X1+X3,d))
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)

## e.lvm2 <- estimate(lvm(Y1[mu:sigma1]~ beta1*X1 + beta2*X2,
##                        Y2[mu:sigma2]~ beta1*X1 + beta2*X2,
##                        Y3[mu:sigma3]~ beta1*X1 + beta2*X2),
##                        data = d)
## e.gls <- gls(value ~ X1 + X2,
##              data = dLred,
##              weight = varIdent(form = ~1|variable),
##              method = "ML")

## test_that("gls equivalent to lvm", {
##     expect_equal(as.double(logLik(e.lvm2)), as.double(logLik(e.gls)))
## })

## ** compare2
test_that("multiple linear regression: df",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    
    name.param <- names(.coef2(e.gls))
    df.gls <- compare2(e.gls, par = name.param, cluster = "Id", bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]

    ## 
    sigma2 <- list(coef(e.lvm)["Y1~~Y1"],
                   coef(e.lvm)["Y2~~Y2"],
                   coef(e.lvm)["Y3~~Y3"])
    X <- list(as.matrix(cbind(1,d[,c("X1")])),
              as.matrix(cbind(1,d[,c("X2")])),
              as.matrix(cbind(1,d[,c("X1","X3")])))
    std.GS <- mapply(X, sigma2, FUN = function(x,y){
        c(sqrt(diag(solve(crossprod(x))*y)),sqrt(2*y^2/n))
    })

    name.coef.lvm <- names(coef(e.lvm))
    expect_equal(df.lvm$std[grep("Y1",name.coef.lvm)], unname(std.GS[[1]]))
    expect_equal(df.lvm$std[grep("Y2",name.coef.lvm)], unname(std.GS[[2]]))
    expect_equal(df.lvm$std[grep("Y3",name.coef.lvm)], unname(std.GS[[3]]))

    ## test value
    df.GS <- lapply(X, function(x){
        c(rep(n,NCOL(x)), n/4)
    })
    expect_equal(df.lvm$df[grep("Y1",name.coef.lvm)], unname(df.GS[[1]]), tol = 1e-7)
    expect_equal(df.lvm$df[grep("Y2",name.coef.lvm)], unname(df.GS[[2]]), tol = 1e-7)
    expect_equal(df.lvm$df[grep("Y3",name.coef.lvm)], unname(df.GS[[3]]), tol = 1e-7)
})

test_that("multiple linear regression: df adjusted",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]
    
    name.param <- names(.coef2(e.gls))
    df.gls <- compare2(e.gls, par = name.param, cluster = "Id", bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]

    ## 
    X <- list(as.matrix(cbind(1,d[,c("X1")])),
              as.matrix(cbind(1,d[,c("X2")])),
              as.matrix(cbind(1,d[,c("X1","X3")])))
    sigma2 <- list(sigma(lm(Y1~X1,data=d))^2,
                   sigma(lm(Y2~X2,data=d))^2,
                   sigma(lm(Y3~X1+X3,data=d))^2)
    std.GS <- mapply(X, sigma2, FUN = function(x,y){
        c(sqrt(diag(solve(crossprod(x))*y)),sqrt(2*y^2/(n-NCOL(x))))
    })

    name.coef.lvm <- names(coef(e.lvm))
    expect_equal(df.lvm$std[grep("Y1",name.coef.lvm)], unname(std.GS[[1]]), tol = 1e-7)
    expect_equal(df.lvm$std[grep("Y2",name.coef.lvm)], unname(std.GS[[2]]), tol = 1e-7)
    expect_equal(df.lvm$std[grep("Y3",name.coef.lvm)], unname(std.GS[[3]]), tol = 1e-7)

    ## test value
    df.GS <- lapply(X, function(x){
        c(rep(n - NCOL(x),NCOL(x)),
          (n - NCOL(x))/4)
    })
    expect_equal(df.lvm$df[grep("Y1",name.coef.lvm)], unname(df.GS[[1]]), tol = 1e-7)
    expect_equal(df.lvm$df[grep("Y2",name.coef.lvm)], unname(df.GS[[2]]), tol = 1e-7)
    expect_equal(df.lvm$df[grep("Y3",name.coef.lvm)], unname(df.GS[[3]]), tol = 1e-7)

})

## * Mixed model: Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~X1+Gender)) 
e.lvm <- estimate(m, d)
## compare2(e.lvm)

e.lmer <- lme4::lmer(value ~ variable + X1 + Gender + (1|Id),
                     data = dLred, REML = FALSE)

e.lme <- nlme::lme(value ~ variable + X1 + Gender, random = ~ 1|Id,
                   data = dLred, method = "ML")

e.gls <- nlme::gls(value ~ variable + X1 + Gender,
                   correlation = corCompSymm(form=~ 1|Id),
                   data = dLred, method = "ML")

## ** clubSandwich - bug
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dLred$Id)
## strange that same type of coef have very different degrees of freedom

## ** compare - ok
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed model: Satterthwaite ",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Satterthwaite")

    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    expect_equal(as.double(GS$df),
                 as.double(df.lvm[1:5,"df"]), tol = 1e-4) ## needed for CRAN
    expect_equal(as.double(GS$tvalue),
                 as.double(df.lvm[1:5,"statistic"]), tol = 1e-8) ## needed for CRAN

    name.param <- names(.coef2(e.lme))
    df.lme <- compare2(e.lme, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    expect_equal(df.lme$statistic, df.lvm$statistic, tol = 1e-5)
    expect_equal(df.lme$df, df.lvm$df, tol = 1e-5)

    name.param <- names(.coef2(e.gls))
    df.gls <- compare2(e.gls, par = name.param, bias.correct = FALSE, as.lava = FALSE)[1:length(name.param),]
    expect_equal(df.gls$statistic[1:5], df.lvm$statistic[1:5], tol = 1e-5)
    expect_equal(df.gls$df[1:5], df.lvm$df[1:5], tol = 1e-5)
})

test_that("mixed model: KR-like correction",{
    ## does not work when running test
    ## GS <- summary(e.lmer, ddf = "Kenward-Roger")$coef[,"df"]
    GS <- lmerTest:::calcSummary(e.lmer, ddf = "Kenward-Roger")

    ## get_Lb_ddf(e.lmer, c(0,1,0,0,0))
    ## get_Lb_ddf(e.lmer, c(0,0,0,1,0))
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]


    name.param <- names(.coef2(e.lme))
    df.lme <- compare2(e.lme, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]
    expect_equal(df.lme$statistic, df.lvm$statistic, tol = 1e-5)
    expect_equal(df.lme$df, df.lvm$df, tol = 1e-5)

    name.param <- names(.coef2(e.gls))
    df.gls <- compare2(e.gls, par = name.param, bias.correct = TRUE, as.lava = FALSE)[1:length(name.param),]
    expect_equal(df.gls$statistic[1:5], df.lvm$statistic[1:5], tol = 1e-5)
    expect_equal(df.gls$df[1:5], df.lvm$df[1:5], tol = 1e-5)
})

## * Mixed model: Unstructured with different variance
m <- lvm(Y1~1*eta,
         Y2~1*eta,
         Y3~1*eta,
         eta~X1+Gender)
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, d)

e.lme <- lme(value ~ variable + X1 + Gender,
             random =~ 1|Id,
             correlation = corSymm(),
             weights = varIdent(form =~ 1|variable),
             data = dLred,
             method = "ML")

e.gls <- gls(value ~ variable + X1 + Gender,
             correlation = corSymm(form=~ 1|Id),
             weights = varIdent(form =~ 1|variable),
             data = dLred,
             method = "ML")

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})

## ** compare
test_that("UN mixed model: df",{
    ## singular information matrix
    ## df.adj.lme <- compare2(e.lme,
    ##                          robust = FALSE, bias.correct = FALSE)
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, par = name.param, bias.correct = FALSE, as.lava = FALSE)

    ## overparametrized model
    ## name.param <- names(.coef2(e.lme))
    ## df.lme <- compare2(e.lme, par = name.param, bias.correct = FALSE, as.lava = FALSE)

    name.param <- names(.coef2(e.gls))
    df.gls <- compare2(e.gls, par = name.param, bias.correct = FALSE, as.lava = FALSE)
})


#----------------------------------------------------------------------
### test-compare2.R ends here


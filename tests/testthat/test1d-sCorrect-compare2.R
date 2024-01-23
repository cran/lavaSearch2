### test-compare2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: jan 23 2024 (11:31) 
##           By: Brice Ozenne
##     Update #: 237
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
## rm(list = ls())
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
d <- lava::sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]

## * linear regression

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
    df.lvm <- compare2(e.lvm, linfct = name.param, ssc = "none", as.lava = FALSE)
        
    ## test value
    n.param <- length(coef(e.lm))
    df.GS <- c(rep(n,n.param), n/4)
    expect_equal(unname(df.lvm$df), df.GS)

    sigma2 <- coef(e.lvm)["Y1~~Y1"]
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/e.lvm$data$n))
    expect_equal(unname(sqrt(diag(df.lvm$vcov))), unname(std.GS))
})

test_that("linear regression: df adjusted",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, linfct = name.param, as.lava = FALSE)

    ## test value
    n.param <- length(coef(e.lvm))-1
    df.GS <- c(rep(n-n.param,n.param), (n-n.param)/4)
    expect_equal(unname(df.lvm$df), df.GS)

    sigma2 <- sigma(e.lm)^2
    iXX <- solve(crossprod(model.matrix(e.lm)))
    std.GS <- c(sqrt(diag(iXX*sigma2)),sqrt(2*sigma2^2/(n-n.param)))
    expect_equal(unname(sqrt(diag(df.lvm$vcov))), unname(std.GS))
})

## * multiple linear regression lvm
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
test_that("multiple linear regression: df adjusted",{
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(estimate2(e.lvm), linfct = name.param, as.lava = FALSE, sep = c("",""))

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
    expect_equal(unname(sqrt(diag(df.lvm$vcov))[c("Y1","Y1~X1","Y1~~Y1")]), unname(std.GS[[1]]), tol = 1e-7)
    expect_equal(unname(sqrt(diag(df.lvm$vcov))[c("Y2","Y2~X2","Y2~~Y2")]), unname(std.GS[[2]]), tol = 1e-7)
    expect_equal(unname(sqrt(diag(df.lvm$vcov))[c("Y3","Y3~X1","Y3~X3","Y3~~Y3")]), unname(std.GS[[3]]), tol = 1e-7)

    ## test value
    df.GS <- lapply(X, function(x){
        c(rep(n - NCOL(x),NCOL(x)),
          (n - NCOL(x))/4)
    })
    expect_equal(unname(df.lvm$df[c("Y1","Y1~X1","Y1~~Y1")]), unname(df.GS[[1]]), tol = 1e-7)
    expect_equal(unname(df.lvm$df[c("Y2","Y2~X2","Y2~~Y2")]), unname(df.GS[[2]]), tol = 1e-7)
    expect_equal(unname(df.lvm$df[c("Y3","Y3~X1","Y3~X3","Y3~~Y3")]), unname(df.GS[[3]]), tol = 1e-7)

})

## * Mixed model: Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~X1+Gender)) 
e.lvm <- estimate(m, d)
## compare2(e.lvm)

## e.lmer <- lmerTest::lmer(value ~ variable + X1 + Gender + (1|Id),
##                           data = dLred, REML = FALSE)
## e2.lmer <- update(e.lmer, REML = TRUE)
e.lme <- nlme::lme(value ~ variable + X1 + Gender, random = ~ 1|Id,
                   data = dLred, method = "ML")

e.gls <- nlme::gls(value ~ variable + X1 + Gender,
                   correlation = corCompSymm(form=~ 1|Id),
                   data = dLred, method = "ML")

## ** clubSandwich - bug
## expect_equal(logLik(e.lmer),logLik(e.lme))
expect_equal(-259.840820322626, as.double(logLik(e.lme)))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dLred$Id)
## strange that same type of coef have very different degrees of freedom

## ** compare 
expect_equal(-259.840820322626,as.double(logLik(e.lvm)))

test_that("mixed model: Satterthwaite ",{

    ## does not work when running test
    ## GS <- summary(e.lmer)$coef[,c("df","t value")]
    GS <- matrix(c(91.83528558, 100.00000017, 100.00000016, 49.99999976, 49.99999977, -1.00221818, 0.63398721, 1.67029848, 10.09702494, 3.08796567), 
                 nrow = 5, 
                 ncol = 2, 
                 dimnames = list(c("(Intercept)", "variableY2", "variableY3", "X1", "GenderFemale"),c("df", "t value")) 
                 ) 

    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, linfct = name.param, ssc = "none", as.lava = FALSE)
    expect_equal(as.double(GS[,"df"]),
                 as.double(df.lvm$df[1:5]), tol = 1e-4) ## needed for CRAN
    expect_equal(as.double(GS[,"t value"]),
                 as.double(summary(df.lvm, test = multcomp::adjusted("none"))$table2[1:5,"statistic"]), tol = 1e-8) ## needed for CRAN

    ## F test
    ## GS <- lmerTest::contestMD(e.lmer, L = diag(1,5,5), rhs = 0, ddf = "Satterthwaite")
    GS <- data.frame("Sum Sq" = c(166.27553944),"Mean Sq" = c(33.25510789),"NumDF" = c(5),"DenDF" = c(76.52281817),"F value" = c(23.33447979),"Pr(>F)" = c(0),
                     check.names = FALSE)
    name.param <- names(coef(e.lvm))    
    df.F <- compare2(e.lvm, linfct = name.param[1:5], ssc = "none", as.lava = FALSE, F.test = TRUE)
    
    expect_equal(GS[["DenDF"]], unname(df.F$global["df"]), tol = 1e-5)
    expect_equal(GS[["F value"]], unname(df.F$global["statistic"]), tol = 1e-8)
})

test_that("mixed model: KR-like correction",{

    ## does not work when running test
    ## GS <- summary(e2.lmer, ddf = "Kenward-Roger")$coef[,c("df","t value")]
    ## get_Lb_ddf(e2.lmer, c(0,1,0,0,0))
    ## get_Lb_ddf(e2.lmer, c(0,0,0,1,0))
    GS <- matrix(c(85.06326979, 98, 98, 47, 47, -0.97751931, 0.62761532, 1.65351114, 9.78942877, 2.99389376), 
                 nrow = 5, 
                 ncol = 2, 
                 dimnames = list(c("(Intercept)", "variableY2", "variableY3", "X1", "GenderFemale"),c("df", "t value")) 
                 ) 
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, linfct = name.param, as.lava = FALSE)

    ## expect_equal(as.double(GS[,"df"]),
    ##              as.double(df.lvm$df[1:5]), tol = 1) ## difference of 1 or 2 df
    expect_equal(as.double(GS[,"t value"]),
                 as.double(summary(df.lvm, test = multcomp::adjusted("none"))$table2[1:5,"statistic"]), tol = 1e-5) ## needed for CRAN

    previous.value <- data.frame("estimate" = c(-0.25588154, 0.15137028, 0.39879913, 1.48076547, 0.92411608, 1.45423356, 0.6594628), 
                                 "se" = c(0.26176621, 0.24118321, 0.24118321, 0.15126168, 0.30866695, 0.20917549, 0.24297288), 
                                 "df" = c(87.2184581, 96.66666667, 96.66666667, 48.33333333, 48.33333333, 24.16666667, 14.29181404), 
                                 "lower" = c(-0.77615186, -0.32733248, -0.07990363, 1.17668766, 0.30361014, 1.02267408, 0.13933428), 
                                 "upper" = c(0.26438879, 0.63007304, 0.87750189, 1.78484328, 1.54462201, 1.88579304, 1.17959131), 
                                 "statistic" = c(-0.97751933, 0.62761531, 1.65351113, 9.78942913, 2.99389387, 6.95221787, 2.71414157), 
                                 "p.value" = c(0.3310157, 0.53173576, 0.10147095, 0, 0.00433193, 3.3e-07, 0.01654271))
    expect_equivalent(previous.value, summary(df.lvm, test = multcomp::adjusted("none"))$table2, tol = 1e-5)
})

### ** compare to SAS
if(FALSE){
    ## setwd("c:/Users/hpl802/AppData/Roaming/R")
    write.table(dLred, file = "mydata.txt", row.names = FALSE)
    ## /* Define path */
    ## %Let NomEtude = %Str(C:\Users\hpl802\AppData\Roaming\R\);

    ## /* Define path to file */
    ## FILENAME Fichier "&NomEtude.%Str(mydata.txt)";

    ## /* Importation of the data */
    ## Data mydata; 
    ## Infile Fichier FirstObs=2 obs=61; /* if no FirstObs : will read all first lines */
    ## input gpr $ animal $ week weight; 
    ## Run;

    ## /* display data */
    ## PROC SGPANEL DATA=mydata;
    ## PANELBY gpr;
    ## SERIES X=week Y=weight / GROUP=animal;
    ## RUN;

    ## /* Fit mixed model */
    ## PROC Mixed DATA = mydata;
    ## Class gpr animal week ;
    ## Model weight = week gpr*week / SOLUTION DDFM=KR;
    ## Repeated week / SUBJECT = animal TYPE = CS R RCORR;
    ## RUN; 
}

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
    name.param <- names(coef(e.lvm))
    df.lvm <- compare2(e.lvm, linfct = name.param, ssc = "none", as.lava = FALSE)

    ##                          estimate       std  statistic       df      p-value
    ## [eta] = 0              -0.2530247 0.2459609 -1.0287194 61.99195 3.076087e-01
    ## [Y2] = 0                0.1513703 0.2248199  0.6732956 50.00000 5.038597e-01
    ## [Y3] = 0                0.3987991 0.2286753  1.7439534 50.00000 8.731285e-02
    ## [eta~X1] = 0            1.4498392 0.1465743  9.8914990 50.00000 2.318146e-13
    ## [eta~GenderFemale] = 0  0.9213738 0.2991017  3.0804696 50.00000 3.355219e-03
    ## [Y1~~Y1] = 0            1.3533853 0.4323206  3.1305133 29.40589 3.923193e-03
    ## [eta~~eta] = 0          0.4391486 0.3092283  1.4201436 21.64808 1.698083e-01
    ## [Y2~~Y2] = 0            1.6200992 0.4392003  3.6887476 13.94845 2.444719e-03
    ## [Y3~~Y3] = 0            1.7889734 0.4646774  3.8499254 13.62954 1.850575e-03
    ## [Y1~~Y2] = 0            0.2231421 0.3296648  0.6768757 24.12389 5.049242e-01
    ## [Y1~~Y3] = 0            0.2638691 0.3376548  0.7814760 23.84905 4.422119e-01
    ## global                         NA        NA 17.0357449 34.39628 7.882273e-11

    previous.value <- data.frame("estimate" = c(-0.25302471, 0.15137028, 0.39879913, 1.4498392, 0.92137382, 1.35338532, 0.43914858, 1.62009916, 1.78897339, 0.22314209, 0.26386911), 
                                 "se" = c(0.24596086, 0.22481994, 0.22867534, 0.14657427, 0.29910174, 0.43232057, 0.30922829, 0.43920032, 0.46467742, 0.32966479, 0.33765478), 
                                 "df" = c(61.99194948, 50, 50, 50, 50, 29.40589311, 21.64808397, 13.94844762, 13.62954237, 24.12389256, 23.8490532), 
                                 "lower" = c(-0.74469474, -0.30019386, -0.0605088, 1.15543612, 0.3206103, 0.46972029, -0.20275661, 0.67778142, 0.78979205, -0.45706773, -0.43324959), 
                                 "upper" = c(0.23864531, 0.60293442, 0.85810706, 1.74424227, 1.52213734, 2.23705035, 1.08105378, 2.5624169, 2.78815472, 0.90335192, 0.9609878), 
                                 "statistic" = c(-1.02871942, 0.6732956, 1.74395338, 9.89149902, 3.08046963, 3.13051333, 1.42014364, 3.68874764, 3.84992537, 0.67687572, 0.78147601), 
                                 "p.value" = c(0.30760871, 0.50385966, 0.08731285, 0, 0.00335522, 0.00392319, 0.16980825, 0.00244472, 0.00185058, 0.50492415, 0.44221192))
    expect_equivalent(previous.value, summary(df.lvm, test = multcomp::adjusted("none"))$table2, tol = 1e-5)
})


#----------------------------------------------------------------------
### test-compare2.R ends here


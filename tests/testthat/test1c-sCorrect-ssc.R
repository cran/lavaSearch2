### test1-sCorrect-ssc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:08) 
## Version: 
## Last-Updated: jan 19 2022 (11:40) 
##           By: Brice Ozenne
##     Update #: 122
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
## rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
library(nlme)
context("sCorrect (small sample correction)")

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
dLred$variable.factor <- as.factor(dLred$variable)

## * linear regression
## ** univariate
e.lm <- lm(Y1~X1+X2, data = d)
e.gls <- gls(Y1~X1+X2, data = d)
e.lvm <- estimate(lvm(Y1~X1+X2), data = d)

test_that("linear regression - residual correction equivalent to REML", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")
    
    ## comapre parameters
    GS <- c(coef(e.lm), sigma(e.lm)^2)
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-6)

    ## compare vcov
    GS <- vcov(e.lm)
    expect_equal(unname(eSSC1.lvm$sCorrect$vcov.param[1:3,1:3]),
                 unname(GS), tol = 1e-6)
        
    ## model based
    GS <- data.frame("estimate" = c(0.15110105, 1.08233491, -0.42993475, 1.82039942), 
                     "se" = c(0.19577922, 0.18717083, 0.20041228, 0.37551973), 
                     "df" = c(47, 47, 47, 11.75), 
                     "lower" = c(-0.24275594, 0.70579577, -0.83311225, 1.00027731), 
                     "upper" = c(0.54495805, 1.45887406, -0.02675726, 2.64052153), 
                     "statistic" = c(0.77179308, 5.78260466, -2.14525157, NA), 
                     "p.value" = c(0.44410038, 5.7e-07, 0.03713074, NA))
    
    expect_equivalent(GS,
                      summary2(eSSC1.lvm)$table2, tol = 1e-6)
    
    ## robust
    GS <- data.frame("estimate" = c(0.15110105, 1.08233491, -0.42993475, 1.82039942), 
                     "se" = c(0.18025537, 0.2334738, 0.15426941, 0.3611838), 
                     "df" = c(47, 47, 47, 11.75), 
                     "lower" = c(-0.21152598, 0.61264622, -0.74028477, 1.03158648), 
                     "upper" = c(0.51372808, 1.55202361, -0.11958474, 2.60921236), 
                     "statistic" = c(0.83826102, 4.63578753, -2.78690869, NA), 
                     "p.value" = c(0.40612711, 2.848e-05, 0.00765275, NA))
    expect_equivalent(GS,
                      summary2(eSSC1.lvm, robust = TRUE)$table2, tol = 1e-6)

})

## test_that("linear regression - Cox correction equivalent to REML", {
##     eSSC2.lvm <- estimate2(e.lvm, ssc = "Cox")

##     ## compare coef
##     GS <- c(coef(e.lm), sigma(e.lm)^2)
##     expect_equal(unname(eSSC2.lm$sCorrect$param),
##                  unname(GS), tol = 1e-6)
##     expect_equal(unname(eSSC2.lvm$sCorrect$param),
##                  unname(GS), tol = 1e-6)
    
##     ## compare vcov
##     GS <- vcov(e.lm)
##     expect_equal(unname(eSSC2.lm$sCorrect$vcov.param[1:3,1:3]),
##                  unname(GS), tol = 1e-6)
##     expect_equal(unname(eSSC2.lvm$sCorrect$vcov.param[1:3,1:3]),
##                  unname(GS), tol = 1e-6)

##     ## compare JJK
##     name.param <- c(names(coef(e.lm)),"sigma2")
##     p <- length(name.param)
##     JJK <- array(0, dim = rep(p,3), dimnames = list(name.param,name.param,name.param))
##     X <- model.matrix(e.lm)

##     JJK[name.param[1:3],name.param[1:3],"sigma2"] <- -crossprod(X)/sigma(e.lm)^4
##     JJK[name.param[1:3],"sigma2",name.param[1:3]] <- -crossprod(X)/sigma(e.lm)^4
##     JJK["sigma2",name.param[1:3],name.param[1:3]] <- crossprod(X)/sigma(e.lm)^4
##     expect_equal(JJK, eSSC2.lm$sCorrect$ssc$JJK, tol = 1e-5)
## })


## ** multiple, no constrain
e.gls <- gls(value ~ -1 + variable + variable:X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "REML")
e.lvm <- estimate(lvm(Y1 ~ X1,
                      Y2 ~ X1,
                      Y3 ~ X1),
                  data = d)

test_that("multiple linear regression - residual correction equivalent to REML", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")

    ## lvm
    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2])*intervals(e.gls)$sigma[2])^2
            )
    
    expect_equal(unname(eSSC1.lvm$sCorrect$param[c("Y1","Y2","Y3","Y1~X1","Y2~X1","Y3~X1","Y1~~Y1","Y2~~Y2","Y3~~Y3")]),
                 unname(GS), tol = 1e-6)

    GS <- vcov(e.gls)
    expect_equal(unname(eSSC1.lvm$sCorrect$vcov.param[1:6,1:6]),
                 unname(GS), tol = 1e-6)

})

## test_that("multiple linear regression - Cox correction equivalent to REML", {
##     eSSC2.lvm <- estimate2(e.lvm, ssc = "Cox")
##     ## eSSC2.gls <- sCorrect(e.gls, ssc = "Cox")
##     ## eSSC2.glsN <- sCorrect(e.gls, ssc = "Cox", derivative = "numeric")

##     ## range(eSSC2.gls$sCorrect$hessian-eSSC2.glsN$sCorrect$hessian)
##     ## eSSC2.gls$sCorrect$information
##     ## eSSC2.glsN$sCorrect$information

##     ## lvm
##     GS <- c(intervals(e.gls)$coef[,2],
##             c(Y1 = 1, intervals(e.gls)$varStruct[,2]^2)*intervals(e.gls)$sigma[2]^2
##             )
    
##     expect_equal(unname(eSSC2.lvm$sCorrect$param[c("Y1","Y2","Y3","Y1~X1","Y2~X1","Y3~X1","Y1~~Y1","Y2~~Y2","Y3~~Y3")]),
##                  unname(GS), tol = 1e-6)
    
##     GS <- vcov(e.gls)
##     expect_equal(unname(eSSC2.lvm$sCorrect$vcov.param[1:6,1:6]),
##                  unname(GS), tol = 1e-6)

##     ## gls
##     ## GS <- c(intervals(e.gls)$coef[,2],
##     ##         sigma2 = as.double(intervals(e.gls)$sigma[2]^2),
##     ##         intervals(e.gls)$varStruct[,2]
##     ##         )
    
## })

## ** multiple, with constrains
e.gls0 <- gls(value ~ variable-1 + X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "ML")
e.gls <- gls(value ~ variable-1 + X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "REML")
e.lvm <- estimate(lvm(Y1[mu1:sigma1]~ beta1*X1,
                      Y2[mu2:sigma2]~ beta1*X1,
                      Y3[mu3:sigma3]~ beta1*X1),
                  data = d)
## logLik(e.gls0)
## logLik(e.gls)
## logLik(e.lvm)

## vcov(e.gls0)[1:4,1:4]/vcov(e.lvm)[1:4,1:4]


test_that("multiple linear regression - residual correction equivalent to REML", {

    expect_equal(information(e.lvm), unname(moments2(e.lvm)$information), tol = 1e-6)
    
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")

    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2]) * intervals(e.gls)$sigma[2])^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-3)
})

## test_that("multiple linear regression - Cox correction equivalent to REML", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

##     GS <- c(intervals(e.gls)$coef[,2],
##             (c(Y1 = 1, intervals(e.gls)$varStruct[,2]) * intervals(e.gls)$sigma[2])^2)

##     ## not precisely the same but better
##     expect_equal(unname(eSSC2.lvm$sCorrect$param),
##                  unname(GS), tol = 1e-3)

## })


## * mixed model
## ** CS
m <- lvm(c(Y1[0:sigma]~1*eta,
           Y2[0:sigma]~1*eta,
           Y3[0:sigma]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   data = dLred, method = "REML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   data = dLred, method = "REML")

test_that("mixed model (CS) - residual correction equivalent to REML", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")
    ## eSSC1.lvm <- sCorrect(update(e.gls, method = "ML"), ssc = "residuals")

    GS <- c(intervals(e.lme)$fixed[,2],
            sigma2 = as.double(intervals(e.lme)$sigma[2])^2,
            tau = intervals(e.lme)$reStruct$Id[,2,]^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-4)
    ## eSSC.lvm$sCorrect$param - GS
})

## test_that("mixed model (CS) - Cox correction equivalent to REML", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
##     ## coef(eSSC2.lvm) - coef(e.lvm)
##     ## GS - coef(e.lvm)
    
##     GS <- c(intervals(e.lme)$fixed[,2],
##             sigma2 = as.double(intervals(e.lme)$sigma[2])^2,
##             tau = intervals(e.lme)$reStruct$Id[,2,]^2)
##     GS2 <- c(intervals(e.gls)$coef[,2],
##             sigma2 = sigma(e.gls)^2,
##             corCoef1 = intervals(e.gls)$corStruct[1,2])

##     ## not precisely the same but better
##     expect_equal(unname(eSSC2.lvm$sCorrect$param),
##                  unname(GS), tol = 1e-4)
##     expect_equal(unname(eSSC2.lme$sCorrect$param),
##                  unname(GS), tol = 1e-4)
##     expect_equal(unname(eSSC2.gls$sCorrect$param),
##                  unname(GS2), tol = 1e-2)
##     ## eSSC.lvm$sCorrect$param - GS
## })

## ** CS with different variances
m <- lvm(c(Y1[0:sigma1]~1*eta,
           Y2[0:sigma2]~1*eta,
           Y3[0:sigma3]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

test_that("mixed model (CS,weight) - residual correction differs from REML", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")

    GS.lme <- c(intervals(e.lme)$fixed[,2],
                Y1 = as.double(intervals(e.lme)$sigma[2])^2,
                tau = intervals(e.lme)$reStruct$Id[,2,]^2,
                intervals(e.lme)$varStruct[,2]^2 * as.double(intervals(e.lme)$sigma[2])^2)

    GS <- c("eta" = 0.37740319, "eta~X1" = 1.32095068, "eta~X2" = -0.02166077, "Y1~~Y1" = 1.07114383, "Y2~~Y2" = 1.50935967, "Y3~~Y3" = 1.78871997, "eta~~eta" = 0.90068904)
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-6)
    ## eSSC.lvm$sCorrect$param - GS
    ## coef(e.lvm) - GS
})

## test_that("mixed model (CS,weight) - Cox correction equivalent to REML", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
##     eSSC2.gls <- sCorrect(e.gls, ssc = "Cox")
##     eSSC2.lme <- sCorrect(e.lme, ssc = "Cox")
##     ## GS - coef(e.lvm)
##     ## coef(eSSC2.lvm) - coef(e.lvm)
    
##     GS <- c(intervals(e.lme)$fixed[,2],
##             Y1 = as.double(intervals(e.lme)$sigma[2])^2,
##             tau = intervals(e.lme)$reStruct$Id[,2,]^2,
##             intervals(e.lme)$varStruct[,2]^2 * as.double(intervals(e.lme)$sigma[2])^2)

##     GS2 <- c(intervals(e.lme)$fixed[,2],
##              c(Y1 = 1,intervals(e.lme)$varStruct[,2])^2 * as.double(intervals(e.lme)$sigma[2])^2,
##              tau = intervals(e.lme)$reStruct$Id[1,2])

##     ## not precisely the same but better
##     expect_equal(unname(eSSC2.lvm$sCorrect$param),
##                  unname(GS), tol = 5e-3)
##     eSSC2.lme$sCorrect$param - GS2
##     ## eSSC.lvm$sCorrect$param - GS
## })


## ** UN
m <- lvm(c(Y1~0+1*eta,
           Y2~0+1*eta,
           Y3~0+1*eta,
           eta~X1+X2))
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, d)

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corSymm(form =~ 1| Id),
                   weight = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

test_that("mixed model (UN) - residual correction differs from REML", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals") 
    ## eSSC1.lvm$sCorrect$param - coef(e.lvm)
    
    gls_sigma2 <- as.double(intervals(e.gls)$sigma[2])^2
    gls_var <- c(Y1 = gls_sigma2, gls_sigma2 * intervals(e.gls)$varStruct[,2]^2)
    gls_tau <- as.double(sqrt(gls_var)["Y2"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[3,2])
    
    ## getVarCov2(eSSC1.lvm)
    ## eSSC1.lvm$sCorrect$param

    GS.gls <- c(intervals(e.gls)$coef[,2],
                Y1 = as.double(gls_var["Y1"]) - gls_tau,
                "eta~~eta" = gls_tau,
                Y2 = as.double(gls_var["Y2"] - gls_tau),
                Y3 = as.double(gls_var["Y3"] - gls_tau),
                "Y1~~Y2" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y2"] * intervals(e.gls)$corStruct[1,2] - gls_tau),
                "Y1~~Y3" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[2,2] - gls_tau)
                )

    GS <- data.frame("estimate" = c(0.40605248, 1.34876061, 0.03171814, 1.36038582, 0.73014581, 1.59213818, 1.8781554, 0.20120537, 0.23243999), 
                     "se" = c(0.16913807, 0.16170109, 0.17314067, 0.48960791, 0.36665137, 0.46517041, 0.50622562, 0.36375786, 0.36974), 
                     "df" = c(48.99131541, 48.99131541, 48.99131541, 28.35663651, 18.85511933, 14.13124925, 13.54461812, 23.82210554, 23.7096419), 
                     "lower" = c(0.06615527, 1.02380864, -0.31622263, 0.35803744, -0.03766363, 0.59531552, 0.78897477, -0.54985076, -0.53116071), 
                     "upper" = c(0.74594969, 1.67371257, 0.37965891, 2.3627342, 1.49795525, 2.58896085, 2.96733602, 0.9522615, 0.99604068), 
                     "statistic" = c(2.40071598, 8.34107309, 0.18319289, NA, NA, NA, NA, 0.55312996, 0.62865795), 
                     "p.value" = c(0.0202054, 0, 0.85540273, NA, NA, NA, NA, 0.5853285, 0.53558303))

    expect_equal(as.double(unlist(model.tables(eSSC1.lvm))),
                 as.double(unlist(GS)), tol = 1e-6)
})

## test_that("mixed model (UN) - Cox correction equivalent to REML", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

##     ## eSSC2.lvm$sCorrect$param - coef(e.lvm)
    
##     gls_sigma2 <- as.double(intervals(e.gls)$sigma[2])^2
##     gls_var <- c(Y1 = gls_sigma2, gls_sigma2 * intervals(e.gls)$varStruct[,2]^2)
##     gls_tau <- as.double(sqrt(gls_var)["Y2"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[3,2])

##     GS <- c(intervals(e.gls)$coef[,2],
##             Y1 = as.double(gls_var["Y1"]) - gls_tau,
##             "eta~~eta" = gls_tau,
##             Y2 = as.double(gls_var["Y2"] - gls_tau),
##             Y3 = as.double(gls_var["Y3"] - gls_tau),
##             "Y1~~Y2" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y2"] * intervals(e.gls)$corStruct[1,2] - gls_tau),
##             "Y1~~Y3" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[2,2] - gls_tau)
##             )

##     ## not precisely the same but better
##     expect_equal(unname(eSSC2.lvm$sCorrect$param),
##                  unname(GS), tol = 5e-3)
##     ## eSSC2.lvm$sCorrect$param - GS
##     ## coef(e.lvm) - GS
## })

## * latent variable model
## ** factor model
m <- lvm(Y1~eta,
         Y2~eta+X2,
         Y3~eta,
         Z1~eta, Z1~~Y1,Z1~~Y2,
         eta~X1+X3)
e.lvm <- estimate(m, d)

## round(coef(estimate(m, sim(m,1e4))),1) ## truth

test_that("factor model - residuals correction", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")
    ## coef(eSSC1.lvm) - coef(e.lvm)

    GS <- data.frame("estimate" = c(0.23990945, 0.28804288, 0.17076884, 0.36590889, 1.1654944, 0.11297345, 0.91569218, 0.52324123, 1.77781303, 0.10836302, 1.19867136, 0.56519387, 1.47128944, 0.28970699, 1.97973709, 0.24899082, 0.35358187), 
                     "se" = c(0.18799645, 0.23072591, 0.2931627, 0.20145141, 0.17621193, 0.10589956, 0.16303206, 0.18016769, 0.22465755, 0.14413976, 0.27267601, 0.18626544, 0.32259805, 0.35063842, 0.40312441, 0.22748203, 0.2547393), 
                     "df" = c(56.44718181, 39.00470949, 18.52136972, 44.292982, 68.78535781, 21.25611374, 12.04379228, 47.13840874, 7.3533839, 17.22586144, 13.16324833, 7.71173505, 12.55326162, 17.06195323, 12.04013252, 22.1993749, 20.8126868), 
                     "lower" = c(-0.13662687, -0.17864252, -0.44390274, -0.04001391, 0.81394165, -0.10709528, 0.56061904, 0.16081869, 1.25171407, -0.19544183, 0.61033243, 0.13285416, 0.77182911, -0.44987084, 1.10172906, -0.2225325, -0.17646796), 
                     "upper" = c(0.61644576, 0.75472828, 0.78544041, 0.77183169, 1.51704715, 0.33304219, 1.27076532, 0.88566378, 2.30391198, 0.41216787, 1.78701029, 0.99753358, 2.17074976, 1.02928481, 2.85774513, 0.72051415, 0.88363171), 
                     "statistic" = c(1.2761382, 1.24842017, 0.58250533, 1.81636304, 6.61416297, 1.06679811, 5.61663881, 2.9041902, 7.91343532, 0.75179128, NA, NA, NA, NA, NA, 1.09455159, 1.38801462), 
                     "p.value" = c(0.20713262, 0.21931911, 0.56725087, 0.07608612, 1e-08, 0.29803039, 0.00011158, 0.00558839, 7.503e-05, 0.46232572, NA, NA, NA, NA, NA, 0.28544743, 0.17981328))

    
    expect_equal(as.double(unlist(model.tables(eSSC1.lvm))),
                 as.double(unlist(GS)),
                 tol = 1e-6)
})

## test_that("factor model - Cox correction", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
##     ## coef(eSSC2.lvm) - coef(e.lvm)
    
##     GS <- c("eta" = 0.23991104, "Y2" = 0.29144175, "Y3" = 0.17832511, "Z1" = 0.36554264, "eta~X1" = 1.16545963, "eta~X3" = 0.11297009, "Y2~eta" = 0.90409861, "Y2~X2" = 0.52324125, "Y3~eta" = 1.75203851, "Z1~eta" = 0.10961229, "Y1~~Y1" = 1.20351764, "eta~~eta" = 0.5603901, "Y2~~Y2" = 1.47834635, "Y3~~Y3" = 0.30540849, "Z1~~Z1" = 1.99182361, "Y1~~Z1" = 0.25025483, "Y2~~Z1" = 0.3555143)
##     expect_equal(coef(eSSC2.lvm),GS, tol = 1e-6)
## })

## ** two factors model (regression)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- lava::estimate(m, d)

test_that("two factors model (correlation) - residuals correction", {
    eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")
    ## summary2(e.lvm, ssc = "residuals")

    GS <- data.frame("estimate" = c(0.1478569, 0.19515562, 0.37384111, 0.39767751, -0.19934296, -0.82545231, 0.36540063, 0.85064787, 0.88015853, 1.29882077, 0.92947602, 1.95754764, 1.22777389, 0.84489487, 2.201649, 1.66547564, 0.8429522, 1.23373245, 0.77442115, 1.40526018, 0.13550962), 
                     "se" = c(0.27670411, 0.21934509, 0.17960657, 0.20040727, 0.24589559, 0.39285267, 0.27933333, 0.17310347, 0.1621181, 0.15302732, 0.26424083, 0.60071584, 0.16064393, 0.3688724, 0.67154771, 0.41296464, 0.30484248, 0.31661237, 0.37777288, 0.33039193, 0.73377475), 
                     "df" = c(24.83880008, 31.95629525, 22.96732922, 48.29904737, 12.20156161, 2.4213353, 7.09306026, 9.4906504, 6.67523697, 47.98657702, 4.53936457, 1.37008364, 47.96401517, 11.36380095, 11.96547193, 13.3573358, 12.98491495, 10.99271212, 6.70939132, 11.93934438, 1.96554411), 
                     "lower" = c(-0.42221348, -0.25165968, 0.00226736, -0.00520386, -0.73412348, -2.26295108, -0.29336504, 0.46212497, 0.49299671, 0.9911365, 0.22895845, -2.17802662, 0.90477135, 0.03617209, 0.73800383, 0.77573987, 0.18430228, 0.53681695, -0.12677344, 0.68499218, -3.07530271), 
                     "upper" = c(0.71792727, 0.64197092, 0.74541486, 0.80055889, 0.33543756, 0.61204646, 1.02416631, 1.23917076, 1.26732036, 1.60650504, 1.62999358, 6.0931219, 1.55077643, 1.65361764, 3.66529418, 2.5552114, 1.50160213, 1.93064796, 1.67561575, 2.12552817, 3.34632195), 
                     "statistic" = c(0.53435021, 0.88971956, 2.08144447, 1.98434679, -0.81068133, -2.10117525, 1.30811684, 4.91410075, 5.42911932, 8.48750911, 3.51753361, 3.25869159, 7.64282757, NA, NA, NA, NA, NA, NA, NA, NA), 
                     "p.value" = c(0.59785007, 0.38026572, 0.04872921, 0.05291443, 0.43307701, 0.1478051, 0.23162872, 0.00071146, 0.00114389, 0, 0.01990062, 0.13317327, 0, NA, NA, NA, NA, NA, NA, NA, NA))

    expect_equal(as.double(unlist(model.tables(eSSC1.lvm))),
                 as.double(unlist(GS)),
                 tol = 1e-6)

})

## test_that("two factors model (correlation) - Cox correction", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

##     GS <- c("eta1" = 0.15334368, "Y2" = 0.19799503, "Y3" = 0.37775788, "eta2" = 0.39767751, "Z2" = -0.18562638, "Z3" = -0.75078602, "eta1~eta2" = 0.35160357, "Y2~eta1" = 0.8409626, "Y3~eta1" = 0.86679841, "Y3~X1" = 1.29882077, "Z2~eta2" = 0.89498429, "Z3~eta2" = 1.76979176, "Z3~X3" = 1.22777389, "Y1~~Y1" = 0.88673054, "eta1~~eta1" = 2.22027337, "Y2~~Y2" = 1.7069357, "Y3~~Y3" = 0.87527537, "Z1~~Z1" = 1.26455176, "eta2~~eta2" = 0.74360185, "Z2~~Z2" = 1.43772669, "Z3~~Z3" = 0.31131063)

##     expect_equal(coef(eSSC2.lvm), GS, tol = 1e-6)
## })

## ## ** two factors model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,eta1~X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,eta2~X2,
           eta1~~eta2))

## e.lvm <- estimate(m, d) ## not done due to lack of convergence
## ## coef(e.lvm)

## test_that("two factors model (correlation) - residuals correction", {
##     eSSC1.lvm <- estimate2(e.lvm, ssc = "residuals")

##     GS <- c("eta1" = 0.1478569, "Y2" = 0.19515562, "Y3" = 0.37384111, "eta2" = 0.39767751, "Z2" = -0.19934296, "Z3" = -0.82545231, "eta1~eta2" = 0.36540063, "Y2~eta1" = 0.85064787, "Y3~eta1" = 0.88015853, "Y3~X1" = 1.29882077, "Z2~eta2" = 0.92947602, "Z3~eta2" = 1.95754764, "Z3~X3" = 1.22777389, "Y1~~Y1" = 0.82498637, "eta1~~eta1" = 2.19604428, "Y2~~Y2" = 1.67543559, "Y3~~Y3" = 0.90484948, "Z1~~Z1" = 1.20086385, "eta2~~eta2" = 0.7944319, "Z2~~Z2" = 1.41920933, "Z3~~Z3" = 0.21553652)
##     expect_equal(coef(eSSC1.lvm),GS, tol = 1e-6)
## })

## test_that("two factors model (correlation) - Cox correction", {
##     eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

##     GS <- c("eta1" = 0.24295941, "Y2" = 0.18388462, "Y3" = 0.32464952, "eta2" = 0.37537296, "Z2" = -0.1851004, "Z3" = -0.77200208, "eta1~X1" = 1.08475547, "Y2~eta1" = 0.88925211, "Y3~eta1" = 1.12263979, "Y3~X1" = 0.82749017, "eta2~X2" = -0.10473803, "Z2~eta2" = 0.89399459, "Z3~eta2" = 1.82282964, "Z3~X3" = 1.21228718, "Y1~~Y1" = 0.98937181, "eta1~~eta1" = 0.9128264, "Y2~~Y2" = 1.59949563, "Y3~~Y3" = 0.75743488, "Z1~~Z1" = 1.29260368, "eta2~~eta2" = 0.72179897, "Z2~~Z2" = 1.46486865, "Z3~~Z3" = 0.22243345, "eta1~~eta2" = 0.14547098)

##     expect_equal(coef(eSSC2.lvm), GS, tol = 1e-6)
##     ## coef(e.lvm) - coef(eSSC2.lvm)
## })

##----------------------------------------------------------------------
### test1-sCorrect-ssc.R ends here

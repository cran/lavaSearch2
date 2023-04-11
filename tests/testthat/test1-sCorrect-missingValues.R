### test1-sCorrect-missingValues.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (13:39) 
## Version: 
## Last-Updated: Jan 17 2022 (23:21) 
##           By: Brice Ozenne
##     Update #: 57
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
context("sCorrect (dealing with missing values)")

## * simulation
n <- 2e1
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
dLred2 <- dL[dL$variable %in% c("Y1","Y2"),]
dLred3 <- dLred2
dLred3[dL$variable == "Y2","Id"] <- paste0("-",dLred3[dL$variable == "Y2","Id"])

## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## using the t test function
e.ttest <- t.test(d$Y1, d$Y2)
e.ttest$parameter

## by hand
sX1 <- var(d$Y1)/n
sX2 <- var(d$Y2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))
df-e.ttest$parameter

## * LVM: factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e.lvm <- estimate(m,d)
e2.lvm <- estimate2(e.lvm)

## ** complete case analysis
missing.Row <- d[1,]
missing.Row[,"Id"] <- -1
missing.Row[,c("Y1","Y2","Y3")] <- NA
## eNA.lvm <- estimate(m, rbind(d,missing.Row), missing = FALSE)
eNA.lvm <- estimate(m, rbind(missing.Row,d,missing.Row), missing = FALSE)

test_that("complete case analysis (factor model)", {
    expect_equal(unname(score2(eNA.lvm, ssc = FALSE, indiv = TRUE)), unname(score(eNA.lvm, indiv = TRUE)))
    ## FIX NA!!!!!!
    
    eNA2.lvm <- estimate2(eNA.lvm)
    
    expect_equal(model.tables(eNA2.lvm), model.tables(e2.lvm))
    ##               estimate        se        df       lower     upper   statistic     p.value
    ## eta1       -0.37387088 0.2962948 18.338440 -0.99554076 0.2477990 -1.26182045 0.222827159
    ## Y2         -0.02252887 0.3297457 15.258888 -0.72432797 0.6792702 -0.06832194 0.946416628
    ## Y3          0.38272845 0.2851559  5.682026 -0.32459821 1.0900551  1.34217267 0.230678148
    ## eta1~X1     0.99599616 0.3134807 18.469394  0.33859531 1.6533970  3.17721651 0.005092730
    ## eta1~X2    -0.04275890 0.2607688  9.859874 -0.62490927 0.5393915 -0.16397243 0.873065310
    ## Y2~eta1     1.05707590 0.2723211  5.481730  0.37515670 1.7389951  3.88172564 0.009722483
    ## Y3~eta1     1.08664682 0.3566308  1.982418 -0.46092704 2.6342207  3.04697992 0.093951240
    ## Y3~X1       0.61495545 0.4296209  3.045199 -0.74088088 1.9707918  1.43139085 0.246430988
    ## Y1~~Y1      0.49861889 0.3398654  5.101312 -0.36983983 1.3670776          NA          NA
    ## eta1~~eta1  1.10299240 0.5611968  4.190452 -0.42760479 2.6335896          NA          NA
    ## Y2~~Y2      1.60214650 0.6249695  4.823447 -0.02223785 3.2265309          NA          NA
    ## Y3~~Y3      0.64437632 0.4273389  4.023035 -0.53943230 1.8281849          NA          NA
})

## ** full information

### *** example of issue with lava
m <- lvm(c(Y1~1*eta1,Y2~1*eta1,Y3~1*eta1+X1,
           eta1~1))

missing.Row <- d[1,]
missing.Row[,"Id"] <- -1
missing.Row[,c("Y1","Y2")] <- NA
missing.Row2 <- d[3,]
missing.Row2[,"Id"] <- -2
missing.Row2[,c("Y1","Y3")] <- NA

dNA.wide <- rbind(missing.Row,d,missing.Row2)
dNA.long <- melt(dNA.wide, measure.vars = c("Y1","Y2","Y3"))
dNA.long$variable <- factor(dNA.long$variable, levels = c("Y1","Y2","Y3"))

test_that("full information (issue lava)", {
    eNA.lvm <- estimate(m, data = dNA.wide, missing = TRUE)
    test <- estimate2(eNA.lvm, ssc = FALSE)

    hessian.GS <- numDeriv::hessian(func = function(x){logLik(eNA.lvm, p = x)},
                                    x = coef(eNA.lvm))
    hessian.info <- information(eNA.lvm, p = coef(eNA.lvm), type = "hessian")
    hessian.lavaSearch2 <- hessian2(test)

    expect_equal(unname(hessian.lavaSearch2), unname(hessian.GS), tol = 1e-6)
    ## expect_equal(unname(hessian.info), unname(hessian.GS), tol = 1e-6) ## fails
})


## *** factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

eNA.lvm <- estimate(m, dNA.wide, missing = TRUE)

test_that("full information (factor model)", {

    hessian.GS <- numDeriv::hessian(func = function(x){logLik(eNA.lvm, p = x)},
                                    x = coef(eNA.lvm))
    expect_equal(hessian.GS, unname(hessian2(eNA.lvm, ssc = FALSE)), tol = 1e-6)

    ## NOT TRUE!!!! issue in lava or on purpose?
    ## expect_equal(unname(information(eNA.lvm)), unname(solve(vcov(eNA.lvm))), tol = 1e-6)
    ## NOT TRUE !!! bug in lava? (see previous ***)
    ## expect_equal(information(eNA.lvm), unname(information2(eNA.lvm, ssc = FALSE)), tol = 1e-6)
    eNA2.lvm <- estimate2(eNA.lvm)
    model.tables(eNA2.lvm)
    ##               estimate        se         df        lower      upper   statistic     p.value
    ## eta1       -0.39593331 0.2854528 19.5265046  -0.99230434  0.2004377 -1.38703618 0.181061270
    ## Y3          0.37704733 0.2887687  8.9341075  -0.27692804  1.0310227  1.30570700 0.224270881
    ## Y2         -0.04308602 0.3263015 16.1994022  -0.73412297  0.6479509 -0.13204356 0.896575991
    ## eta1~X1     1.03194106 0.2972545 17.5733098   0.40634382  1.6575383  3.47157405 0.002802916
    ## eta1~X2    -0.01552331 0.2462704  9.4842975  -0.56831932  0.5372727 -0.06303358 0.951048335
    ## Y3~eta1     1.07977655 0.3656814  1.0724740  -2.88655281  5.0461059  2.95277987 0.194108317
    ## Y3~X1       0.62193065 0.4431914  2.2225543  -1.11318213  2.3570434  1.40330027 0.283921347
    ## eta1~~eta1  1.02409942 0.5259689  0.5684192 -44.51400617 46.5622050          NA          NA
    ## Y3~~Y3      0.62918891 0.4126684  4.8415855  -0.44213402  1.7005118          NA          NA
    ## Y2~eta1     1.07161261 0.2676486  8.1125205   0.45590084  1.6873244  4.00380494 0.003818464
    ## Y1~~Y1      0.51334322 0.3337466  2.9800835  -0.55281455  1.5795010          NA          NA
    ## Y2~~Y2      1.54266511 0.6031601  5.4456513   0.02959063  3.0557396          NA          NA
    expect_equal(summary(eNA2.lvm)$table2$df,
                 c(19.52650464, 8.93410748, 16.19940217, 17.57330981, 9.4842975, 1.07247398, 2.2225543, 0.56841916, 4.84158554, 8.1125205, 2.9800835, 5.44565132),
                 tol = 1e-6)
})

##----------------------------------------------------------------------
### test1-sCorrect-missingValues.R ends here

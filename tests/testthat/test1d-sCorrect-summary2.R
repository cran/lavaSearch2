### test1-sCorrect-summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  4 2018 (13:29) 
## Version: 
## Last-Updated: jan 17 2022 (16:51) 
##           By: Brice Ozenne
##     Update #: 67
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

    printDF <- function(object, bias.correct){
        colDF <- summary2(object, bias.correct = bias.correct)$coef[,"df",drop=FALSE]
        n.coef <- NROW(colDF)
        vec.end <- c(rep(",",n.coef-1),")")
        vec.start <- c("c(", rep("",n.coef-1))        
        df <- data.frame(paste0(vec.start,"\"",rownames(colDF),"\""),
                         "=",
                         paste0(colDF[,1],vec.end))
        names(df) <- NULL
        print(df, row.names = FALSE)
    }

}

lava.options(symbols = c("~","~~"))
library(nlme)
context("sCorrect (dVcov-SatterthwaiteCorrection)")

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
## ** model fit
e.lvm <- estimate(lvm(Y1~X1+X2+Gender), data = d)

## ** test df
test_that("linear regression: Satterthwaite (df)", {
    df <- c("Y1" = 50,
            "Y1~X1" = 50,
            "Y1~X2" = 50,
            "Y1~GenderFemale" = 50,            
            "Y1~~Y1" = 12.5)
    expect_equal(as.double(df),
                 summary2(e.lvm, ssc = "none")$coef$df)
})

test_that("linear regression: Satterthwaite + SSC (df)", {
    df <- c("Y1" =   46,
            "Y1~X1" =   46,
            "Y1~X2" =   46,
            "Y1~GenderFemale" =   46,
            "Y1~~Y1" = 11.5)
    expect_equal(as.double(df),
                 summary2(e.lvm, ssc = "residuals")$coef$df)
})

## ** robust standard error
test_that("linear regression: robust SE", {
    ## printDF(e.lvm, bias.correct = TRUE)
    lava.options(df.robust = 1)
    eS0 <- summary2(e.lvm, robust = TRUE, df = "satterthwaite")$coef
    
    df <- c("Y1" =   46,
            "Y1~X1" =   46,
            "Y1~X2" =   46,
            "Y1~GenderFemale" =   46,
            "Y1~~Y1" = 11.5)
    expect_equal(as.double(df),
                 eS0$df, tol = 1e-2)
    
    eS1 <- summary2(e.lvm, robust = TRUE, df = "satterthwaite")$coef
    eS2 <- summary2(e.lvm, robust = TRUE, df = "satterthwaite", cluster = 1:n)$coef
    expect_equal(eS1,eS2)
   
})

## * linear regression with constrains 
## ** model fit
e.lvm <- estimate(lvm(Y1[0:2]~X1+1*X2), data = d)

e.lvm2 <- estimate(lvm(Y1~beta*X1+beta*X2), d)


## ** test df
test_that("linear regression with constrains: Satterthwaite (df)", {
    expect_equal(summary2(e.lvm)$coef$df,c(Inf)) ## Inf since the variance coefficient is known
    ## printDF(e.lvm2, bias.correct = FALSE)
    df <- c("Y1~X1" =   50,
            "Y1" =   50,
            "Y1~~Y1" = 12.5)
    expect_equal(summary2(e.lvm2, ssc = "none")$coef$df,
                 as.double(df))
})

test_that("linear regression with constrains: Satterthwaite + SSC (df)", {
    expect_equal(summary2(e.lvm)$coef$df,c(Inf)) ## Inf since the variance coefficient is known
    ## printDF(e.lvm2, bias.correct = TRUE)
    df <- c("Y1~X1" =   48,
            "Y1" =   48,
            "Y1~~Y1" = 12)
    expect_equal(summary2(e.lvm2)$coef$df,
                 as.double(df))
})

## * multiple linear regression 
## ** model fit
ls.lm <- list(lm(Y1~X1,d),lm(Y2~X2,d),lm(Y3~X1+X3,d))
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)

## ** test df
test_that("multiple linear regression: Satterthwaite (df)", {
    ## printDF(e.lvm, bias.correct = FALSE)
    df <- c("Y1~X1" =   50,
            "Y2~X2" =   50,
            "Y3~X1" =   50,
            "Y3~X3" =   50,
            "Y1~~Y1" = 12.5,
            "Y2~~Y2" = 12.5,
            "Y3~~Y3" = 12.5,
            "Y1" =   50,
            "Y2" =   50,
            "Y3" =   50)
    expect_equal(summary2(e.lvm, ssc = "none")$coef[names(df),"df"],
                 as.double(df)) ## 
    
})

test_that("multiple linear regression: Satterthwaite + SSC (df)", {
    ## printDF(e.lvm, bias.correct = TRUE)
    df <- c("Y1~X1" =    48,
            "Y2~X2" =    48,
            "Y3~X1" =    47,
            "Y3~X3" =    47,
            "Y1~~Y1" =    12,
            "Y2~~Y2" =    12,
            "Y3~~Y3" = 11.75,
            "Y1" =    48,
            "Y2" =    48,
            "Y3" =    47)
    expect_equal(summary2(e.lvm)$coef[names(df),"df"],
                 as.double(df)) ## 
    
})

## * multiple linear regression with constrains 
## ** model fit
e.lvm <- estimate(lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2), data = d)

## ** test df
test_that("multiple linear regression with constrains: Satterthwaite (df)", {
    ## printDF(e.lvm, bias.correct = FALSE)
    df <- c("Y1~X1" =   50,
            "Y1~X2" =   NA,
            "Y2~X1" =   NA,
            "Y2~X3" =   NA,
            "Y3~X2" =   50,
            "Y1~~Y1" = 12.5,
            "Y2~~Y2" = 12.5,
            "Y3~~Y3" = 12.5,
            "Y1" =   50,
            "Y2" =   50,
            "Y3" =   50)
    expect_equal(summary2(e.lvm, ssc = "none")$coef[names(df),"df"],
                 as.double(df)) ## 
    
})

test_that("multiple linear regression with constrains: Satterthwaite + SSC (df)", {
    ## printDF(e.lvm, bias.correct = TRUE)
    df <- c("Y1~X1" =    48,
            "Y1~X2" =    NA,
            "Y2~X1" =    NA,
            "Y2~X3" =    NA,
            "Y3~X2" =    48,
            "Y1~~Y1" =    12,
            "Y2~~Y2" = 12.25,
            "Y3~~Y3" =    12,
            "Y1" =    48,
            "Y2" =    49,
            "Y3" =    48)
    expect_equal(summary2(e.lvm)$coef[names(df),"df"],
                 as.double(df)) ## 
    
})

## * multiple linear regression with covariance links 
## ** model fit
e.lvm <- estimate(lvm(Y1~X1+X2,Y2~X3+X1,Y3~X2,Y1~~Y2),d)

## ** test df
test_that("multiple linear regression with covariance: Satterthwaite (df)", {
    ## printDF(e.lvm, bias.correct = FALSE)
    df <- c("Y1~X1" = 50.0023249929247,
            "Y1~X2" = 50.0557533452502,
            "Y2~X1" = 50.1412333709522,
            "Y2~X3" = 50.0557533452502,
            "Y3~X2" =               50,
            "Y1~~Y1" =             12.5,
            "Y1~~Y2" = 14.4382586892588,
            "Y2~~Y2" =             12.5,
            "Y3~~Y3" =             12.5,
            "Y1" = 51.0449669789772,
            "Y2" = 50.0000667169911,
            "Y3" =               50)
    expect_equal(summary2(e.lvm, ssc = "none")$coef[names(df),"df"],
                 as.double(df), tol = 1e-6) ## 
    
})

test_that("multiple linear regression with covariance: Satterthwaite +SSC (df)", {
    ## printDF(e.lvm, bias.correct = TRUE)
    df <- c("Y1~X1" = 47.0021511585814,
            "Y1~X2" = 47.0515840309539,
            "Y2~X1" = 47.1306527469107,
            "Y2~X3" = 47.0515840309539,
            "Y3~X2" =               48,
            "Y1~~Y1" =            11.75,
            "Y1~~Y2" = 13.6588307493286,
            "Y2~~Y2" =            11.75,
            "Y3~~Y3" =               12,
            "Y1" = 47.9656506208306,
            "Y2" = 47.0000617288762,
            "Y3" =               48)
    expect_equal(summary2(e.lvm)$coef[names(df),"df"],
                 as.double(df), tol = 1e-6) ## 
    
})



## * two factors model (regression)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- lava::estimate(m, d)

test_that("Two factor model: Satterthwaite +SSC (df)", {
    GS <- model.tables(estimate2(e.lvm))
    expect_equal(summary2(e.lvm)$coef[rownames(GS),"df"],
                 GS$df)
})
######################################################################
### test1-sCorrect-summary2.R ends here

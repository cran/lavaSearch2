### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: feb  4 2018 (13:55) 
##           By: Brice Ozenne
##     Update #: 142
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## source("c:/Users/hpl802/Documents/GitHub/lavaSearch2/tests/testthat/test1-iid2-lava.R")

## * header
if(FALSE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

library(clubSandwich)
lava.options(symbols = c("~","~~"))

context("iid2-lava")

n <- 5e1

## * linear model
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

set.seed(10)
m <- lvm(formula.lvm)
transform(m,Id~Y) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(m,n)

e.lm <- lm(formula.lvm,data=d)
e.lvm <- estimate(lvm(formula.lvm),data=d)
name.coef <- names(coef(e.lm))

## ** iid2 matches iid
test_that("iid2 matches iid", {
    e.iid2.lm <- iid2(e.lm, bias.correct = FALSE) ## dim(e.iid2.lm)
    GS1 <- iid(e.lm) ## dim(GS1) ## head(GS1)
    attr(GS1, "bread") <- NULL
    expect_equal(e.iid2.lm[,name.coef], GS1)

    e1.iid2.lvm <- iid2(e.lvm, bias.correct = FALSE)
    expect_equal(unname(e1.iid2.lvm), unname(e.iid2.lm))
    e.GS <- iid(e.lvm)
    attr(e.GS, "bread") <- NULL

    ## NOTE: iid in lava uses numerical derivative to compute the information matrix
    ## this is why there is not a perfect matching between iid2.lvm and iid.lvm
    I <- numDeriv::jacobian(function(p){
        score(e.lvm, p = p, indiv = FALSE)
    }, pars(e.lvm), method = lava.options()$Dmethod)
    score.GS <- score(e.lvm, indiv = TRUE)

    expect_equivalent(e1.iid2.lvm, score.GS %*% vcov(e.lvm))
    expect_equivalent(e.GS, - score.GS %*% solve(I))
})


## ** iid2 lvm matches iid2 lm
test_that("iid2 lvm matches iid2 lm", {
    for(iAdj in c(FALSE,TRUE)){ # iAdj <- 1
        e.iid2.lm <- iid2(e.lm, bias.correct = iAdj)
        e0.iid2.lvm <- iid2(e.lvm, bias.correct = iAdj)
        
        expect_equivalent(e.iid2.lm, e0.iid2.lvm)       
    }
})

## ** iid2 matches clubSandwich
test_that("iid2.lm/iid2.lvm matches clubSandwich", {
    eHC2.iid2.lm <- iid2(e.lm, bias.correct = TRUE, as.clubSandwich = 2)
    ## Aleardy checked before
    ## eHC2.iid2.lvm <- iid2(e.lvm, bias.correct = TRUE, as.clubSandwich = 2)
    ## expect_equal(unname(eHC2.iid2.lm),
    ##              unname(eHC2.iid2.lvm))

    VsandwichHC2.lm <- crossprod(eHC2.iid2.lm)

    V.GS <- clubSandwich::vcovCR(e.lm, type = "CR2", cluster = d$Id)
    expect_equal(as.matrix(V.GS),
                 VsandwichHC2.lm[name.coef,name.coef])
})




#----------------------------------------------------------------------
### test-iid2.R ends here

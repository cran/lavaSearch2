### demoGLS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  6 2018 (16:49) 
## Version: 
## Last-Updated: feb  6 2018 (17:04) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

# devtools::install_github("bozenne/lavaSearch2")
library(lavaSearch2)
library(nlme)

## * simulate data
mSim <- lvm(c(Y1,Y2,Y3) ~ eta, eta ~ 0.1 * Treatment)
categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
transform(mSim, Id ~ eta) <- function(x){1:NROW(x)}
latent(mSim) <- ~eta

set.seed(10)
dfW.data <- sim(mSim, n = 2e1, latent = FALSE)

## * data management: from wide to long format
dfL.data <- reshape2::melt(dfW.data, id.vars = c("Treatment","Id"))
dfL.data <- dfL.data[order(dfL.data$Id,dfL.data$variable),]

## * fit gls model
e.gls <- gls(value ~ variable + Treatment,
             correlation = corSymm(form =~ 1|Id),
             weights = varIdent(form =~ 1|variable),
             data = dfL.data)

## Wald test according to nlme:::gls
summary(e.gls)$tTable

## Wald test according to lavaSearch2
summary2(e.gls, bias.correct = TRUE)$tTable

## NOTE: each time you call summary2 it will re-compute the small sample correction
e.gls2 <- e.gls
sCorrect(e.gls2) <- TRUE # equivalent to bias.correct = TRUE
summary(e.gls2)$tTable # fast

## F test
compare2(e.gls2, par = c("variableY2 = 0", "variableY3 = 0"))
# or via a contrast matrix
resC <- createContrast(e.gls2, par = c("variableY2 = 0", "variableY3 = 0"),
                       add.variance = TRUE)
compare2(e.gls2, contrast = resC$contrast, null = resC$null)

## * fit lme model
e.lme <- lme(value ~ variable + Treatment,
             random = ~ 1|Id,
             correlation = corSymm(form =~ 1|Id),
             weights = varIdent(form =~ 1|variable),
             data = dfL.data)

## same model as gls
logLik(e.lme)-logLik(e.gls)

## Wald test according to nlme:::gls
summary(e.lme)$tTable ## better small sample correction for treatment (df = 18 instead of 38)

## Wald test according to lavaSearch2
summary2(e.lme, bias.correct = TRUE)$tTable

## NOTE: each time you call summary2 it will re-compute the small sample correction
e.lme2 <- e.lme
sCorrect(e.lme2) <- TRUE # compute the small sample correction
summary(e.lme2)$tTable # fast


##----------------------------------------------------------------------
### demoGLS.R ends here

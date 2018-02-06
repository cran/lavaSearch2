### demoGLS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  6 2018 (16:49) 
## Version: 
## Last-Updated: feb  6 2018 (16:49) 
##           By: Brice Ozenne
##     Update #: 2
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

#### simulate data ####
mSim <- lvm(c(Y1,Y2,Y3) ~ eta, eta ~ 0.1 * Treatment)
categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
transform(mSim, Id ~ eta) <- function(x){1:NROW(x)}
latent(mSim) <- ~eta

set.seed(10)
dfW.data <- sim(mSim, n = 2e1, latent = FALSE)

#### from wide to long format
dfL.data <- reshape2::melt(dfW.data, id.vars = c("Treatment","Id"))
dfL.data <- dfL.data[order(dfL.data$Id,dfL.data$variable),]

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
sCorrect(e.gls2) <- TRUE # compute the small sample correction
summary(e.gls2)$tTable # fast


##----------------------------------------------------------------------
### demoGLS.R ends here

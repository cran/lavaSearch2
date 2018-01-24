### Lava_modelsearchLR.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (17:58) 
## Version: 
## last-updated: jan 18 2018 (14:19) 
##           By: Brice Ozenne
##     Update #: 83
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - modelsearchLR
#' @title Testing the Relevance of Additional Links Using a Likelihood Ratio Test
#' @description Testing the Relevance of Additional Links Using a Likelihood Ratio Test.
#' @name modelsearchLR
#' 
#' @return a list
#'
#' @keywords internal

## * modelsearchLR
#' @rdname modelsearchLR
modelsearchLR <- function (x, data, restricted, link, directive, 
                           update.FCT, update.args,
                           method.p.adjust, display.warnings, trace){

    ## ** initialisation
    n.link <- length(link)
    df.test <- data.frame("link" = link,
                          "statistic" = as.numeric(rep(NA,n.link)),
                          "p.value" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value" = as.numeric(rep(NA,n.link)),
                          "convergence" = as.numeric(rep(NA,n.link)),
                          "coefBeta" = as.numeric(rep(NA,n.link)),
                          "corrected.level" = as.numeric(rep(NA,n.link)),
                          "quantile" = as.numeric(rep(NA,n.link)),
                          stringsAsFactors = FALSE)

    best.test <- -Inf
    best.model <- NULL
        
    if(trace > 0){pb <- utils::txtProgressBar(max = n.link, style = 3) }

    for (iterI in 1:n.link) { # iterI <- 1
        newfit <- update.FCT(x, args = update.args,
                             restricted = restricted[iterI,], directive = directive[iterI])

        if(class(newfit) != "try-error" && !is.na(stats::logLik(newfit))){ 

            if(newfit$opt$convergence == 0 ){ # test whether the model has correctly converged
                newCoef.tempo <- stats::coef(newfit)[setdiff(names(coef(newfit)),names(coef(x)))]
                df.test[iterI, "coefBeta"] <- newCoef.tempo
                if(class(newfit) == "lvmfit"){
                    compareT <- lava::compare(x,newfit)
                    df.test[iterI, "statistic"] <- compareT$statistic[[1]]
                    df.test[iterI, "p.value"] <- compareT$p.value[[1]]
                }else{
                    compareT <- stats::anova(x, newfit)
                    df.test[iterI, "statistic"] <- compareT$F[2]
                    df.test[iterI, "p.value"] <- compareT$`Pr(>F)`[2]
                }
                df.test[iterI, "convergence"] <- 0
            }else{
                df.test[iterI, "convergence"] <- 1             
            }
 
            if(!is.na(df.test[iterI,"statistic"]) && df.test[iterI,"statistic"]>best.test){
                best.test <- df.test[iterI,"statistic"]
                best.model <- newfit
            }
        }    
    
        if(trace > 0){ utils::setTxtProgressBar(pb, value = iterI) }    
    }
    if(trace > 0){  close(pb) }
    df.test$adjusted.p.value <- stats::p.adjust(df.test$p.value, method = method.p.adjust)
    
    ## ** export 
    return(list(df.test = df.test,
                best.test = best.test,
                best.model = best.model))
}

#----------------------------------------------------------------------
### Lava_modelsearchLR.R ends here

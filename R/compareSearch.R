### compareSearch.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (11:57) 
## Version: 
## last-updated: mar 22 2018 (17:26) 
##           By: Brice Ozenne
##     Update #: 304
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - compareSearch
#' @title Compare Methods to Identify Missing Local Dependencies in a LVM
#' @description Compare methods to identify missing local dependencies in a LVM.
#' @name compareSearch
#' 
#' @param object a \code{lvm} model.
#' @param alpha [numeric 0-1] the significance cutoff for the p-values.
#' When the p-value is below, the corresponding link will be added to the model
#' and the search will continue. Otherwise the search will stop.
#' @param statistic [character] statistic used to perform the test.
#' Can the likelihood ratio test (\code{"LR"}),
#' the score (\code{"score"}),
#' or the max statistic (\code{"max"}).
#' @param method.p.adjust [character] the method used to adjust the p.values for multiple comparisons.
#' Can be any method that is valid for the \code{stats::p.adjust} function (e.g. \code{"fdr"}).
#' Ignored when using the max statistic.
#' @param trace [logical] should the execution of the function be traced?
#' @param ... arguments passed to \code{\link{modelsearch2}}.
#'
#' @details This function calls the \code{\link{modelsearch2}} function
#' to find the local dependencies.
#' 
#' @return A list containing:
#' \itemize{
#' \item newlink: a list containing for each \code{statistic}-\code{method.p.adjust} the local dependencies.
#' \item table.coef: a \code{data.frame} object containing for each \code{statistic}-\code{method.p.adjust} the estimated coefficients.
#' \item ls.search: a list containing for each \code{statistic}-\code{method.p.adjust} a \code{modelsearch2} object.
#' }
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#'
#' \dontrun{
#' res <- compareSearch(e.lvm, statistic = c("score","Wald"),
#'                      method.p.adjust = c("holm","fdr","max"))
#' res
#' }
#' @concept modelsearch

## * function - compareSearch
#' @rdname compareSearch
#' @export
compareSearch <- function(object, alpha = 0.05,
                          method.p.adjust, statistic,
                          trace = 1, ...){
    
    p.value <- link <- NULL
    
    ### ** normalize arguments
    method.p.adjust <- sapply(method.p.adjust, match.arg, choices = lava.options()$search.p.adjust, several.ok = TRUE)
    statistic <-  sapply(statistic, match.arg, choices = lava.options()$search.statistic, several.ok = TRUE)

    ### ** fit all models using unadjusted p.values
    ls.search <- list()
    if("score" %in% statistic){
        if(trace){
            cat("modelsearch with the score statistic")
        }
        ls.search$score <- try(modelsearch2(object, statistic = "score", method.p.adjust = "none",
                                            trace = trace-1, ...), silent = TRUE)

        if("try-error" %in% class(ls.search$score)){
            statistic <- setdiff(statistic,"score")
            if(trace){
                cat(" - error \n")
            }
        }else if(trace){
            cat(" - done \n")
        }
    }
    if("LR" %in% statistic){
        if(trace){
            cat("modelsearch with the log likelihood ratio statistic")
        }
        ls.search$LR <- try(modelsearch2(object, statistic = "LR", method.p.adjust = "none",
                                         trace = trace-1, ...), silent = TRUE)
        if("try-error" %in% class(ls.search$LR)){
            statistic <- setdiff(statistic,"LR")
            if(trace){
                cat(" - error \n")
            }
        }else if(trace){
            cat(" - done \n")            
        }
    }
    if("Wald" %in% statistic){
        if(trace){
            cat("modelsearch with the robust Wald statistic")
        }
        if(any(c("fastmax","max") %in% method.p.adjust)){
            if("fastmax" %in% method.p.adjust){
                ls.search$Wald <- try(modelsearch2(object, statistic = "Wald", method.p.adjust = "fastmax",
                                                   trace = trace-1, ...), silent = TRUE)
               method.p.adjust[method.p.adjust == "fastmax"] <- "max"
            }else{
                ls.search$Wald <- try(modelsearch2(object, statistic = "Wald", method.p.adjust = "max",
                                                   trace = trace-1, ...), silent = TRUE)
            }
            

            ## *** check whether further steps are needed for other type 1 error adjustment
            ## i.e. it could be that after max adjustement the p.value is >0.05 but if we don't adjust we should continue
            if("try-error" %in% class(ls.search$Wald) == FALSE){
                currentStep <- nStep(ls.search$Wald)
                vec.tempo <- getStep(ls.search$Wald, step = currentStep, slot = "sequenceTest")
                maxStep <- list(...)$nStep
                vec.p.adjust <- sapply(setdiff(method.p.adjust,c("fastmax","max")), function(iAdj){
                    min(stats::p.adjust(vec.tempo$p.value, method = iAdj))
                })
                if(is.null(maxStep)){maxStep <- Inf}

                ## *** perform additional search
                if(any(vec.p.adjust < alpha) && (vec.tempo$selected==FALSE) && (currentStep<maxStep) ){ # continue the modelsearch

                    ## add the link of the last test to the model (avoid to repeat step)
                    model.tempo <- getStep(ls.search$Wald, step=nStep(ls.search$Wald), slot = "sequenceModel")
                    link.tempo <- getStep(ls.search$Wald, step=nStep(ls.search$Wald), slot = "sequenceTest")
                    newLink.tempo <- link.tempo[which.min(p.value),link]

                    ls.args <- lapply(model.tempo$call[-(1:2)], evalInParentEnv)
                    restricted.tempo <- unlist(initVarLink(newLink.tempo))
                    directive.tempo <- length(grep(lava.options()$symbols[2],newLink.tempo,fixed=TRUE))==0
                
                    if("lvmfit" %in% class(object)){
                        model.tempo2 <- .updateModelLink.lvm(model.tempo, args = ls.args,
                                                             restricted = restricted.tempo,
                                                             directive = directive.tempo)
                    }else{
                        model.tempo2 <- .updateModelLink.default(model.tempo, args = ls.args,
                                                                 restricted = restricted.tempo,
                                                                 directive = directive.tempo)
                    }

                    dots <- list(...)
                    dots$nStep <- maxStep-currentStep
                    if("link" %in% names(dots)){
                        dots$link <- setdiff(dots$link, union(getNewLink(ls.search$Wald, step = 1:currentStep),newLink.tempo))
                    }
                    otherSearch <- try(do.call("modelsearch2", c(list(x = model.tempo2,
                                                                      statistic = "Wald",
                                                                      method.p.adjust = "none",
                                                                      trace = trace-1), dots)
                                               ), silent = TRUE)

                    if("try-error" %in% class(otherSearch) == FALSE){
                        ls.search$Wald <- merge(ls.search$Wald, otherSearch)
                    }
                }
             }
            
        }else{            
            ls.search$Wald <- try(modelsearch2(object, statistic = "Wald", method.p.adjust = "none", 
                                               trace = FALSE, ...), silent = TRUE)            
        }
        if("try-error" %in% class(ls.search$Wald)){
            statistic <- setdiff(statistic,"Wald")
            if(trace){
                cat(" - error \n")
            }
        }else if(trace){
            cat(" - done \n")
        }
    }

    ### ** check only errors
    if(length(statistic)==0){
        return(list(newlinks = NULL,
                    table.coef = NULL,
                    ls.search = ls.search)
               )
    }

### ** Adjust p.values
    ls.searchAll <- list()    
    for(iStatistic in statistic){ # iStatistic <- statistic[1]
        ##  print(iStatistic)
        for(iAdjust in method.p.adjust){ # iAdjust <- method.p.adjust[1]
            ## print(iAdjust)
            if(iAdjust == "max" && iStatistic != "Wald"){next}
            list.tempo <- list(.adjustModelSearch(ls.search[[iStatistic]],
                                                  model0 = object,
                                                  method.p.adjust = iAdjust,
                                                  alpha  = alpha))
            names(list.tempo) <- paste0(iStatistic,"-",iAdjust)
            ls.searchAll <- c(ls.searchAll,list.tempo)
        }
            
    }

### ** Merge results
    ## newlinks
    ls.newlinks <- lapply(ls.searchAll,getNewLink)

    ## value of all links
    name.alllinks <- unique(unlist(lapply(ls.searchAll, function(x){
        names(coef(getStep(x, step = nStep(x), slot = "sequenceModel")))
    })))
    name.search <- names(ls.searchAll)
    table.alllinks <- matrix(NA, nrow = length(name.alllinks), ncol = length(ls.searchAll)+1,
                             dimnames = list(name.alllinks, c("base",name.search)))
    table.alllinks[names(stats::coef(object)),"base"] <- coef(object)
    for(iSearch in name.search){ # iSearch <- name.search[2]
        M.tempo <- getStep(ls.searchAll[[iSearch]], step = nStep(ls.searchAll[[iSearch]]), slot = "sequenceModel")
        table.alllinks[names(stats::coef(M.tempo)),iSearch] <- coef(M.tempo)
    }
    
    ### ** export
    return(list(newlinks = ls.newlinks,
                table.coef = table.alllinks,
                ls.search = ls.searchAll))
    
}

## * adjustModelSearch
.adjustModelSearch <- function(object, model0,  method.p.adjust, alpha){

    n.Test <- length(object$sequenceTest)
    
    ## ** adjust p.value
    seqP.value <- rep(NA, n.Test)
    for(iTest in 1:n.Test){
        if(method.p.adjust!="max"){            
            object$sequenceTest[[iTest]]$adjusted.p.value <- stats::p.adjust(object$sequenceTest[[iTest]]$p.value,
                                                                             method = method.p.adjust)
        }
        seqP.value[iTest] <- min(object$sequenceTest[[iTest]]$adjusted.p.value)
    }

    ## ** stop search when necessary
    index.keepTest <- union(1, which(seqP.value<alpha)+1)
    index.keepTest <- index.keepTest[sapply(index.keepTest, function(x){
        all(1:x %in% index.keepTest) # remove non consecutive steps
    })]
    index.keepTest <- index.keepTest[index.keepTest<=nStep(object)] # remove extra step due to early stop

    seqP.value <- seqP.value[index.keepTest]
    object$sequenceTest <- object$sequenceTest[index.keepTest]
    object$sequenceModel <- object$sequenceModel[index.keepTest]
    if(method.p.adjust=="max"){ # max never activated
        object$sequenceQuantile <- object$sequenceQuantile[index.keepTest]
        object$sequenceIID <- object$sequenceIID[index.keepTest]
        object$sequenceSigma <- object$sequenceSigma[index.keepTest]
    }

    ## ** update final model
    index.finalModel <- utils::tail(which(seqP.value<alpha),1)
    if( length(index.finalModel) == 0 ){
        object$sequenceModel[[1]] <- model0
    }else{
        object$sequenceModel[[length(object$sequenceModel)]] <- object$sequenceModel[[index.finalModel]]
    }
    
    ## ** update adjustment
    object$method.p.adjust <- method.p.adjust
    
    ## ** export
    return(object)
}

#----------------------------------------------------------------------
### compareSearch.R ends here

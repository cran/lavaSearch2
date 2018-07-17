### summary.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (10:46) 
## Version: 
## last-updated: maj  2 2018 (10:41) 
##           By: Brice Ozenne
##     Update #: 73
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * method summary.modelsearch2
#' @title summary Method for modelsearch2 Objects
#' @description summary method for modelsearch2 objects.
#'
#' @param object output of the \code{modelsearch2} function.
#' @param print should the summary be printed in the terminal.
#' @param ... [internal] only used by the generic method.
#' 
#' @method summary modelsearch2
#' @export
summary.modelsearch2 <- function(object, print = TRUE, ...){

    p.value <- NULL # [:for CRAN check] subset
    
    ## ** extract data from object
    xx <- object$sequenceTest
    df.seqTest <- do.call(rbind,lapply(xx, function(step){
        if("convergence" %in% names(step)){
            step$noConvergence <- sum(step$convergence!=0)
            step$convergence <- sum(step$convergence==0)
        }
        indexMax <- which.max(abs(step$statistic))
        return(step[indexMax,,drop = FALSE])
    }))
    n.step <- NROW(df.seqTest)
    n.selected <- sum(df.seqTest$selected)

    keep.cols <- c("link","nTests","noConvergence","statistic","adjusted.p.value")    
    if(!is.na(object$method.p.adjust) && object$method.p.adjust == "max"){
        keep.cols <- c(keep.cols,"quantile")
    }

    if("df" %in% names(df.seqTest)){
        df.seqTest$nTests.adj <- 0.05/(2*(1-stats::pt(df.seqTest$quantile,
                                                      df = df.seqTest$df)))
    }else{
        df.seqTest$nTests.adj <- 0.05/(2*(1-stats::pnorm(df.seqTest$quantile)))
    }
    if(object$method.p.adjust == "fastmax"){
        df.seqTest[p.value==0, c("p.value", "adjusted.p.value")] <- NA
    }
    
    ## ** output
    out <- list(output = list(), data = df.seqTest)
    statistic <- switch(object$statistic,
                        "Wald" = "Wald",
                        "score" = "score",
                        "LR" = "likelihood ratio",
                        "NA" = "NA")
    if(statistic=="Wald"){
        addOn <- switch(object$typeSD,
                        "information"="",
                        "robust"="robust ",
                        "jackknife"="jackknife ")
        statistic <- paste0(addOn, statistic)
    }
    
    out$output$message.pre <- paste0("Sequential search for local dependence using the ",statistic," statistic \n")
    if(n.selected==0){
        out$output$message.pre <- c(out$output$message.pre,
                                    "The variable selection procedure did not retain any variable \n")
    }else{
        out$output$message.pre <- c(out$output$message.pre,
                                    paste0("The variable selection procedure retained ",n.selected," variable",
                                           if(n.selected>1){"s"},":\n")
                                    )     
    }
     
    out$output$table <- df.seqTest[,keep.cols,drop=FALSE]
    out$output$message.post <- paste0("confidence level: ",1-object$alpha," (two sided, adjustement: ",object$method.p.adjust,")\n")  

    ## ** display
    if(print){
        cat(out$output$message.pre)
        print(out$output$table)
        cat(out$output$message.post)        
    }
    
    ## ** export
    return(invisible(out))
}



#----------------------------------------------------------------------
### summary.modelsearch2.R ends here

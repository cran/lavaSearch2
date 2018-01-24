### print.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 25 2017 (10:18) 
## Version: 
## last-updated: sep 18 2017 (11:16) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.modelsearch2
#' @method print modelsearch2
#' @export
print.modelsearch2 <- function(x, ...){

    out <- summary(x, display = TRUE, ...)
    return(invisible(out))
}

#----------------------------------------------------------------------
### print.modelsearch2.R ends here

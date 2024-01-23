### p.adjust2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 19 2019 (11:28) 
## Version: 
## Last-Updated: Jan 11 2022 (17:38) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


p.adjust2 <- function(p, method, vcov.param = NULL){

    traditional.method <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    new.method <- c("AB1","AB2")
    all.method <- c(traditional.method,new.method)
    
    if(length(method)!=1){
        stop("Argument \'method\' must have length 1 \n")
    }
    
    if(method %in% traditional.method){
        out <- stats::p.adjust(p = p, method = method)
    }else{
        if(method %in% new.method == FALSE){
            stop("Argument \'method\' must be one of \"",paste(all.method,collapse = "\" \""),"\"\n")
        }
        
        if(method %in% c("AB1","AB2") && is.null(vcov.param)){
            stop("Argument \'vcov.param\' must not be NULL when argument \'method\' is \"AB1\" or \"AB2\" \n")
        }
        if(method %in% c("AB1","AB2")){
            
            if(is.null(names(p))){
                stop("Argument \'p\' must be named \n")
            }
            if(!all(names(p) %in% colnames(vcov.param))){
                stop("The column names of argument \'vcov.param\' must match the names of argument \'p\' \n")
            }

            ## compute average correlation
            M.rho <- stats::cov2cor(vcov.param)
            diag(M.rho) <- NA
            r <- colMeans(abs(M.rho[,names(p),drop=FALSE]), na.rm = TRUE)
            k <- length(p)
            out <- switch(method,
                          "AB1" = pmin(1, p*(k-(k-1)*sqrt(abs(r)))),
                          "AB2" = pmin(1, p*k^(1-sqrt(abs(r))))
                          )
            attr(out,"r") <- r
                   
            
        }
    }

    names(out) <- names(p)

    return(out)
}

######################################################################
### p.adjust2.R ends here

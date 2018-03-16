### Lava_modelsearchMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (18:32) 
## Version: 
## last-updated: mar 12 2018 (17:18) 
##           By: Brice Ozenne
##     Update #: 701
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - modelsearchMax
#' @title Testing the Relevance of Additional Links Using the Max Statistic
#' @description Testing the Relevance of Additional Links Using the Max Statistic.
#' 
#' @name modelsearchMax
#'
#' @return A \code{lvmfit} object.
#'
#' @seealso \code{link{modelsearch2}}
#'
#' @concept modelsearch
#' @keywords internal
#'


## * Function - modelsearchMax
## ' @rdname modelsearchMax
modelsearchMax <- function(x, restricted, link, directive, packages,
                           update.FCT, update.args, iid.FCT,                           
                           alpha, method.p.adjust, method.max, n.sim = 1e3, 
                           iid.previous = NULL, quantile.previous = NULL, 
                           export.iid, trace, ncpus, init.cpus){

    ## WARNING: do not put link as NULL for data.table since it is used as an argument by the function
    
    ### ** initialisation
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    n.link <- NROW(restricted)
    nObs <- NROW(update.args$data)

    best.test <- -Inf
    best.model <- NULL    
    iid.link <- NULL
    convergence <- rep(NA,n.link)

    typeSD <- attr(iid.FCT, "typeSD")
    df <- attr(iid.FCT, "df")
    bias.correct <- attr(iid.FCT, "bias.correct")
    
    ### ** wraper
    warper <- function(iterI){ # iterI <- 2

        out <- list(df = data.frame(statistic = as.numeric(NA),
                                    df = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    adjusted.p.value = as.numeric(NA),
                                    convergence = as.numeric(NA),
                                    coefBeta = as.numeric(NA),
                                    stringsAsFactors = FALSE),
                    iid = NULL)
        ## *** fit new model
        newfit <- update.FCT(x, args = update.args,
                             restricted = restricted[iterI,], directive = directive[iterI])
        out$df[1, "convergence"] <- newfit$opt$convergence
        
        ## *** extract influence function        
        if(class(newfit) != "try-error"){ # test whether the model was estimated
            if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged

                ## extract coefficient
                new.coef <- stats::coef(newfit)
                if(link[iterI] %in% names(new.coef) == FALSE){
                    stop("Coefficient ",link[iterI]," not found \n",
                         "Possible coefficients: ",paste0(names(new.coef), collapse = " "),"\n")
                }

                out$df[1, "coefBeta"] <- new.coef[link[iterI]]
                ## extract degree of freedom and standard error
                if(df || bias.correct){
                    sCorrect(newfit, df = df, score = TRUE) <- bias.correct
                    out$iid <- iid2(newfit, robust = (typeSD != "information") )[,link[iterI],drop=FALSE]
                    if(df){
                        e.df <- compare2(newfit, par = link[iterI], as.lava = FALSE)
                        out$df[1, "df"] <- e.df[1, "df"]
                    }
                    sd.coef <- sqrt(sum(out$iid^2, na.rm = TRUE))                    
                }else{
                    out$iid <- sqrt(nObs)*iid.FCT(newfit)[,link[iterI],drop=FALSE]
                    if(typeSD == "information"){
                        ## NOTE: it assumes that whenever df is FALSE then bias.correct is FALSE
                        sd.coef <- sqrt(stats::vcov(newfit)[link[iterI],link[iterI]])
                        out$iid <- out$iid * sd.coef / sd(out$iid, na.rm = TRUE)
                    }else{
                        sd.coef <- sqrt(sum(out$iid^2, na.rm = TRUE))
                    }
                        
                }

                ## compute test statistic
                out$df[1, "statistic"] <- abs(out$df$coefBeta/sd.coef) ## keep .SD for clarity
            }
        }
        return(out)
    }
    
### ** get influence function
    if(trace>0){
        cat("gather influence functions \n")
    }
            
    if(ncpus>1){

        FCTcombine <- function(res1,res2){
            res <- list(df = rbind(res1$df,res2$df),
                        iid = cbind(res1$iid,res2$iid))
            return(res)
        }

        if(init.cpus){
            test.package <- try(requireNamespace("doParallel"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'doParallel\' \n",
                     "This package is necessary when argument \'ncpus\' is greater than 1 \n")
            }
            test.package <- try(requireNamespace("foreach"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'foreach\' \n",
                     "This package is necessary when argument \'ncpus\' is greater than 1 \n")
            
            }
            cl <- parallel::makeCluster(ncpus)
            doParallel::registerDoParallel(cl)
        }
    
        if(trace > 0){
            test.package <- try(requireNamespace("tcltk"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'tcltk\' \n",
                     "This package is necessary when argument \'trace\' is TRUE \n")
            }
            parallel::clusterExport(cl, varlist = "trace")
        }

        vec.packages <- c("lavaSearch2", packages)
        i <- NULL # [:for CRAN check] foreach
        res <- foreach::`%dopar%`(
                            foreach::foreach(i = 1:n.link, .packages =  vec.packages,
                                             # .export = c("ls.LVMargs"),
                                             .combine = FCTcombine),
                            {
                                if(trace){
                                    if(!exists("pb")){
                                        pb <- tcltk::tkProgressBar("modelsearchMax:", min=1, max=n.link)
                                    }
                                    tcltk::setTkProgressBar(pb, i)
                                }
                                return(warper(i))
                            })

        if(init.cpus){
            parallel::stopCluster(cl)
        }
        
    }else{
        if(trace>0){
            test.package <- try(requireNamespace("pbapply"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'pbapply\' \n",
                     "This package is necessary when argument \'trace\' is TRUE \n")
            }                
            resApply <- pbapply::pblapply(1:n.link, warper)
            
        }else{
            resApply <- lapply(1:n.link, warper)
        }
        res <- list(df = do.call(rbind, lapply(resApply,"[[","df")),
                    iid = do.call(cbind,lapply(resApply,"[[","iid")))
        
    }
    df.test <- cbind(link = link, res$df)    
    iid.link <- res$iid

    if(all(df.test$convergence!=0)){
        stop("none of the extended model has converged \n",
             "the additional links may be misspecified \n")
    }
    
    ### ** p.value
    indexCV <- which(df.test$convergence==0)
    df.test[indexCV, "p.value"] <- as.numeric(NA)
    if(df){
        df.test[indexCV, "p.value"] <- 2*(1-stats::pt(abs(df.test[indexCV,"statistic"]),
                                                      df = df.test[indexCV,"df"]))
    }else{
        df.test[indexCV, "p.value"] <- 2*(1-pnorm(abs(df.test[indexCV,"statistic"])))
    }

    ### ** adjust p.value
    if(method.p.adjust == "fastmax"){

        adj.tempo <- stats::p.adjust(df.test$p.value, method = "bonferroni")
        if(any(adj.tempo<0.05)){
            df.test$p.value <- as.numeric(df.test$p.value != min(df.test$p.value))
            method.p.adjust <- "bonferroni"
        }else{
            method.p.adjust <- "max"
        }
        
    }
    
    if(method.p.adjust == "max"){
        nameN0 <- df.test[indexCV, "link"]
        statisticN0 <- setNames(subset(df.test, subset = convergence==0, select = "statistic", drop = TRUE),
                                nameN0)
        if(df){
            dfN0.all <- subset(df.test, subset = convergence==0, select = "df", drop = TRUE)
            dfN0 <- round(stats::median(dfN0.all))
        }else {
            dfN0 <- NULL
        }
        resQmax <- calcDistMaxIntegral(statistic = statisticN0[1], iid = iid.link[,1,drop=FALSE], df = dfN0,
                                       iid.previous = iid.previous, quantile.previous = quantile.previous, 
                                       alpha = alpha, ncpus = ncpus, init.cpus = FALSE, trace = trace)
        
        if(method.max=="integration"){
            resQmax <- calcDistMaxIntegral(statistic = statisticN0,
                                           iid = iid.link,
                                           df = dfN0,
                                           iid.previous = iid.previous,
                                           quantile.previous = quantile.previous, 
                                           alpha = alpha,
                                           ncpus = ncpus,
                                           init.cpus = FALSE,
                                           trace = trace)
            resQmax$p.adjust
        }else{
            method.boot <- switch(method.max,
                                  "boot-naive" = "naive",
                                  "boot-residual" = "residual",
                                  "boot-wild" = "wild")
            
            resQmax <- calcDistMaxBootstrap(statistic = statisticN0, iid = iid.link, method = method.boot, n.sim = n.sim,
                                            iid.previous = iid.previous, quantile.previous = quantile.previous, 
                                            alpha = alpha, ncpus = ncpus, init.cpus = FALSE, trace = trace)
        }
        df.test[indexCV, "corrected.level"] <- resQmax$correctedLevel
        df.test[indexCV, "adjusted.p.value"] <- resQmax$p.adjust
        df.test[indexCV, "quantile"] <- resQmax$z
        Sigma <- resQmax$Sigma
        rownames(Sigma) <- df.test[indexCV, "link"]
        colnames(Sigma) <- df.test[indexCV, "link"]
        
        if(init.cpus){
            parallel::stopCluster(cl)
        }
        
    }else{
        df.test[indexCV, "corrected.level"] <- as.numeric(NA)
        df.test[indexCV, "adjusted.p.value"] <- stats::p.adjust(df.test[indexCV, "p.value"],
                                                                method = method.p.adjust)
        df.test[indexCV, "quantile"] <- as.numeric(NA)
        Sigma <- NULL        
    }    

    ### ** export
    out <- list(df.test = df.test,
                iid = if(export.iid){iid.link}else{NULL},
                Sigma = Sigma)
    return(out)
}

## ----------------------------------------------------------------------
## Lava_modelsearchMax.R ends here

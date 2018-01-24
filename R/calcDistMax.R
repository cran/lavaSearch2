### calcDistMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 21 2017 (16:44) 
## Version: 
## last-updated: jan 22 2018 (11:40) 
##           By: Brice Ozenne
##     Update #: 410
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation
#' @title Adjust the p.values Using the Quantiles of the Max Statistic
#' @description Adjust the p.values using the quantiles of the max statistic.
#' @name calcDistMax
#'
#' @param statistic the observed statistic relative to the coefficients to test.
#' @param iid zero-mean iid decomposition of the observed coefficients used to compute the statistic.
#' @param iid.previous zero-mean iid decomposition of the previous step to condition on.
#' @param quantile.compute should the critical quantile be computed.
#' @param quantile.previous critical threshold of the previous step to condition on.
#' If not \code{NULL} the values should correspond the variable in to the first column(s) of the argument iid.
#' @param df the degree of freedom for the t statistic.
#' @param method the method used to compute the p.values. Can be \code{"integration"}, \code{"boot-wild"}, or \code{"boot-norm"}.
#' See the detail section.
#' @param alpha the significance threshold for retaining a new link
#' @param ncpus the number of cpu to use for parellel computations
#' @param initCpus should the cpus be initialized.
#' @param n.sim the total number of simulations.
#' @param n.repMax the maximum number of rejection when using "\code{"boot-wild"} or \code{"boot-norm"}.
#' @param trace should the execution of the function be traced.
#' 
#' @examples 
#' library(mvtnorm)
#'
#' set.seed(10)
#' n <- 100
#' p <- 4
#' link <- letters[1:p]
#' n.sim <- 1e3 # number of bootstrap simulations 
#'
#' #### test - not conditional ####
#' X.iid <- rmvnorm(n, mean = rep(0,p), sigma = diag(1,p))
#' colnames(X.iid) <- link
#' statistic <- setNames(1:p,link)
#'
#' 
#' r1 <- calcDistMaxIntegral(statistic = statistic, iid = X.iid, 
#'             trace = FALSE, alpha = 0.05, df = 1e6) 
#' 
#' r2 <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             method = "naive",
#'             trace = FALSE, alpha = 0.05, n.sim = n.sim)
#'
#' r3 <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             method = "residual",
#'             trace = FALSE, alpha = 0.05, n.sim = n.sim)
#'
#' r4 <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             method = "wild",
#'             trace = FALSE, alpha = 0.05, initCpus = TRUE, n.sim = n.sim)
#' 
#' rbind(integration = c(r1$p.adjust, quantile = r1$z),
#'       bootNaive    = c(r2$p.adjust, quantile = r2$z),
#'       bootResidual = c(r3$p.adjust, quantile = r3$z),
#'       bootWild    = c(r4$p.adjust, quantile = r4$z))
#'
#' #### test - conditional ####
#' \dontrun{
#' Z.iid <- rmvnorm(n, mean = rep(0,p+1), sigma = diag(1,p+1))
#' seqQuantile <- qmvnorm(p = 0.95, delta = rep(0,p+1), sigma = diag(1,p+1), 
#'                     tail = "both.tails")$quantile
#' 
#' r1c <- calcDistMaxIntegral(statistic = statistic, iid = X.iid,
#'             iid.previous = Z.iid, quantile.previous =  seqQuantile, 
#'             trace = FALSE, alpha = 0.05, df = NULL)
#' 
#' r2c <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             iid.previous = Z.iid, quantile.previous =  seqQuantile, method = "naive",
#'             trace = FALSE, alpha = 0.05, n.sim = n.sim)
#' 
#' r3c <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             iid.previous = Z.iid, quantile.previous =  seqQuantile, method = "residual",
#'             trace = FALSE, alpha = 0.05, n.sim = n.sim)
#'
#' r4c <- calcDistMaxBootstrap(statistic = statistic, iid = X.iid,
#'             iid.previous = Z.iid, quantile.previous =  seqQuantile, method = "wild",
#'             trace = FALSE, alpha = 0.05, n.sim = n.sim)
#'
#' rbind(integration = c(r1c$p.adjust, quantile = r1c$z),
#'       bootNaive    = c(r2c$p.adjust, quantile = r2c$z),
#'       bootResidual = c(r3c$p.adjust, quantile = r3c$z),
#'       bootWild    = c(r4c$p.adjust, quantile = r4c$z))
#' }
#' 


## * calcDistMaxIntegral
#' @rdname calcDistMax
#' @export
calcDistMaxIntegral <- function(statistic, iid, df, 
                                iid.previous = NULL, quantile.previous = NULL,
                                quantile.compute = lava.options()$search.calc.quantile.int,
                                alpha, ncpus = 1, initCpus = TRUE, trace){

    ## ** normalize arguments
    p.iid <- NCOL(iid)
    n <- NROW(iid)
    conditional <- length(quantile.previous)
    if(length(quantile.previous)>1){
        stop("Can only condition on one previous step \n")
    }
    if(is.null(df)){
        distribution.statistic <- "gaussian"
    }else{
        distribution.statistic <- "student"
    }
    
    iid.all <- cbind(iid,iid.previous)
    index.new <- 1:NCOL(iid)
    index.previous <- setdiff(1:NCOL(iid.all),index.new)
    p.iid.all <- NCOL(iid.all)

    ## ** Compute the correlation matrix between the test statistics
    # center to be under the null
    # scale since we want the distribution of the Wald statistic (i.e. statistic with unit variance)
    iid.statistic <- scale(iid.all, center = TRUE, scale = TRUE)
    Sigma.statistic <- stats::cov(iid.statistic, use = "pairwise.complete.obs")
    out <- list(p.adjust = NULL, z = NULL, Sigma = Sigma.statistic[index.new,index.new,drop=FALSE])

    ## ** Definition of the functions used to compute the quantiles
    warperQ <- function(alpha){
        .calcQmaxIntegration(alpha = alpha, p = p.iid,
                             Sigma = Sigma.statistic[index.new,index.new,drop=FALSE],
                             df = df, distribution = distribution.statistic)
    }
    warperP <- function(index){
        .calcPmaxIntegration(statistic = statistic[index], p = p.iid,
                             Sigma = Sigma.statistic[index.new,index.new,drop=FALSE],
                             df = df, distribution = distribution.statistic)
    }
        
    ## ** correction for conditioning on the previous steps
    if(conditional==TRUE){
        out$correctedLevel <- calcType1postSelection(1-alpha, quantile.previous =  quantile.previous,
                                                      mu = rep(0,p.iid.all), Sigma = Sigma.statistic,
                                                      distribution =  distribution.statistic,
                                                      df = df)
        alpha <- 1-out$correctedLevel
    }else{
        out$correctedLevel <- NA
    }

    ## ** Computation
    if(quantile.compute){       
        out$z <- warperQ(alpha)
    }else{
        out$z <- NA
    }
    
    if(trace > 0){ cat("Computation of multivariate student probabilities to adjust the p.values: ") }
    if(ncpus > 1){
        ## *** parallel computations
        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doParallel::registerDoParallel(cl)
        }

        if(trace > 0){
            pb.max <- length(index.new)
            parallel::clusterExport(cl, "trace")
        }

        value <- NULL # [:for CRAN check] foreach
        out$p.adjust <- foreach::`%dopar%`(
                                     foreach::foreach(value = index.new,
                                                      .packages = c("tmvtnorm","mvtnorm"),
                                                      .export = c(".calcPmaxIntegration"),
                                                      .combine = "c"),
                                     {
                                         if(trace){
                                             if(!exists("pb")){
                                                 pb <- tcltk::tkProgressBar("calcDistMaxIntegral:", min=1, max=pb.max)
                                             }
                                             tcltk::setTkProgressBar(pb, value)
                                         }
                                         return(warperP(value))
                                     })

        if(initCpus){
            parallel::stopCluster(cl)
        }
            
    }else{
        ## *** sequential computations
        if(trace>0){
            requireNamespace("pbapply")
            out$p.adjust <- pbapply::pbsapply(index.new, warperP)
        }else{
            out$p.adjust <- sapply(index.new, warperP)
        }
                        
    }
    out$p.adjust <- stats::setNames(out$p.adjust, names(statistic))
        
    if(trace > 0){ cat("done \n") }

    ## ** export
    return(out)
}

## * calcDistMaxBootstrap
#' @rdname calcDistMax
#' @export
calcDistMaxBootstrap <- function(statistic, iid, iid.previous = NULL, quantile.previous = NULL,
                                 method, alpha, ncpus = 1, initCpus = TRUE, n.sim, trace, n.repMax = 100){

    ## ** normalize arguments
    n <- NROW(iid)
    conditional <- length(quantile.previous)>0
    if(length(quantile.previous)>1){
        stop("Can only condition on one previous step \n")
    }

    iid.all <- cbind(iid,iid.previous)
    index.new <- 1:NCOL(iid)
    index.previous <- setdiff(1:NCOL(iid.all),index.new)

    ## ** Function used for the simulations
    warperBoot <- .bootMaxDist
    
    ## ** Compute the correlation matrix between the test statistics
    # center to be under the null
    # scale since we want the distribution of the Wald statistic (i.e. statistic with unit variance)
    iid.statistic <- scale(iid.all, center = TRUE, scale = TRUE)
    Sigma.statistic <- stats::cov(iid.statistic, use = "pairwise.complete.obs")

    ## ** Computation
    if(trace > 0){ cat("Bootsrap simulations to get the 95% quantile of the max statistic: ") }

    if(ncpus>1){
        n.simCpus <- rep(round(n.sim/ncpus),ncpus)
        n.simCpus[1] <- n.sim-sum(n.simCpus[-1])

        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doParallel::registerDoParallel(cl)
        }
  
        i <- NULL # [:for CRAN check] foreach
        distMax <- foreach::`%dopar%`(
                                foreach::foreach(i = 1:ncpus, .packages =  c("MASS"),
                                                 .export = "calcDistMax",
                                                 .combine = "c"),{
                                                     replicate(n.simCpus[i],
                                                               warperBoot(iid = iid.all, sigma = Sigma.statistic,
                                                                          n = n, method = method,
                                                                          index.new = index.new, index.previous = index.previous,
                                                                          quantile.previous = quantile.previous, n.repMax = n.repMax))
                                                 })

        if(initCpus){
            parallel::stopCluster(cl)
        }
        
    }else{

       if(trace>0){
            requireNamespace("pbapply")
            distMax <- pbapply::pbsapply(1:n.sim, warperBoot, method = method,
                                         iid = iid.all, sigma = Sigma.statistic, n = n,                                         
                                         index.new = index.new, index.previous = index.previous,
                                         quantile.previous = quantile.previous, n.repMax = n.repMax)
       }else{
           distMax <- sapply(1:n.sim, warperBoot, method = method,
                             iid = iid.all, sigma = Sigma.statistic, n = n,
                             index.new = index.new, index.previous = index.previous,
                             quantile.previous = quantile.previous, n.repMax = n.repMax)
       }
        
    }
     
    if(trace > 0){ cat("done \n") }

    ## ** export
    out <- list()
    out$z <- stats::quantile(distMax, probs = 1-alpha, na.rm = TRUE)
    out$p.adjust <- sapply(abs(statistic), function(x){mean(distMax>x,na.rm=TRUE)})
    out$Sigma <- Sigma.statistic
    out$correctedLevel <- NA
    return(out)
}
    
## * .calcQmaxIntegration: numerical integration to compute the critical threshold
.calcQmaxIntegration <- function(alpha, p, Sigma, df, distribution){

    if(distribution == "gaussian"){
        if(p==1){
            q.alpha <- stats::qnorm(1-alpha, mean = 0, sd = 1)
        }else{
            q.alpha <- mvtnorm::qmvnorm(1-alpha,
                                        mean = rep(0,p),
                                        corr = Sigma,
                                        tail = "both.tails")$quantile
        }
    }else if(distribution == "student"){
        if(p==1){
            q.alpha <- stats::qt(1-alpha, df = df)
        }else{
            q.alpha <- mvtnorm::qmvt(1-alpha,
                                     delta = rep(0,p),
                                     corr = Sigma,
                                     df = df,
                                     tail = "both.tails")$quantile
        }
    }

    return(q.alpha)
}
    
## * .calcPmaxIntegration_firstStep: numerical integration to compute the p.values
.calcPmaxIntegration <- function(statistic, p, Sigma, df, distribution){
    value <- abs(statistic)
    if(!is.na(value)){
        if(distribution == "gaussian"){
            if(p==1){
                p <- stats::pnorm(value, mean = 0, sd = Sigma)-stats::pnorm(-value, mean = 0, sd = Sigma)
            }else{                
                p <- mvtnorm::pmvnorm(lower = -value, upper = value,
                                      mean = rep(0, p), corr = Sigma)
            }
        }else if(distribution == "student"){
            if(p==1){
                p <- stats::pt(value, df = df)-stats::pt(-value, df = df)
            }else{
                p <- mvtnorm::pmvt(lower = -value, upper = value,
                                   delta = rep(0, p), corr = Sigma, df = df)
            }
        }
        return(1-p)
    }else{
        return(NA)
    }   
}

## * .bootMaxDist: bootstrap simulation
.bootMaxDist <- function(iid, sigma, n, method,
                         index.new, index.previous, quantile.previous, n.repMax,
                         ...){

    iRep <- 0
    cv <- FALSE

    while(iRep < n.repMax && cv == FALSE){

        ## ** resample to obtain a new influence function
        if(method == "naive"){
            iid.sim <- iid[sample.int(n, replace = TRUE),]
        }else if(method == "residual"){
            iid.sim <- MASS::mvrnorm(n,rep(0,NCOL(sigma)),sigma)                    
        }else if(method == "wild"){
            e <- stats::rnorm(n,mean=0,sd=1)
            iid.sim <- sapply(1:NCOL(sigma),function(x){e*iid[,x]})        
        }
        if(!is.null(quantile.previous)){
            iid.previous <- iid.sim[,index.previous]
            test.previous <- apply(iid.previous,2,function(x){sqrt(n)*mean(x)/stats::sd(x)})
            max.previous <- max(abs(test.previous))        
            if(max.previous<quantile.previous){
                iRep <- iRep + 1
            }else{
                iid.sim <- iid.sim[,index.new]
                cv <- TRUE
            }
        }else{
            cv <- TRUE
        }
    }
    
    ## ** compute the bootstrap test statistic
    if(cv){
        Test <- apply(iid.sim,2,function(x){sqrt(n)*mean(x)/stats::sd(x)})
    }else{
        Test <- NA
    }
    return(max(abs(Test)))
}


#----------------------------------------------------------------------
### calcDistMax.R ends here

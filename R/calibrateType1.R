### calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (10:23) 
## Version: 
## Last-Updated: maj 28 2018 (23:47) 
##           By: Brice Ozenne
##     Update #: 518
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Documentation - calibrateType1
##' @title Simulation Study Assessing Bias and Type 1 Error
##' @description Perform a simulation study over one or several sample size
##' to assess the bias of the estimate
##' and the type 1 error of the Wald test and robust Wald test
##' @name calibrateType1
##' 
##' @param object a \code{lvm} object defining the model to be fitted.
##' @param null [character vector] names of the coefficient whose value will be tested against 0. 
##' @param n [integer vector, >0] sample size(s) considered in the simulation study.
##' @param n.rep [integer, >0] number of simulations per sample size.
##' @param cluster  [integer vector] the grouping variable relative to which the observations are iid.
##' Will be passed to \code{lava::estimate}.
##' @param generative.object [lvm] object defining the statistical model generating the data.
##' @param generative.coef [name numeric vector] values for the parameters of the generative model.
##' Can also be \code{NULL}: in such a case the coefficients are set to default values decided by lava (usually 0 or 1).
##' @param true.coef [name numeric vector] expected values for the parameters of the fitted model.
##' @param n.true [integer, >0] sample size at which the estimated coefficients will be a reliable approximation of the true coefficients.
##' @param round.true [integer, >0] the number of decimal places to be used for the true value of the coefficients. No rounding is done if \code{NULL}.
##' @param bootstrap [logical] should bootstrap resampling be performed?
##' @param n.bootstrap [integer, >0] the number of bootstrap sample to be used for each bootstrap.
##' @param checkType1 [logical] returns an error if the coefficients associated to the null hypotheses do not equal 0.
##' @param checkType2 [logical] returns an error if the coefficients associated to the null hypotheses equal 0.
##' @param dir.save [character] path to the directory were the results should be exported.
##' Can also be \code{NULL}: in such a case the results are not exported.
##' @param F.test [logical] should a multivariate Wald test be perform testing simultaneously all the null hypotheses?
##' @param label.file [character] element to include in the file name.
##' @param seed [integer, >0] seed value that will be set at the beginning of the simulation to enable eproducibility of the results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param cpus [integer >0] the number of processors to use.
##' If greater than 1, the simulations are performed in parallel. 
##' @param trace [integer] should the execution of the function be trace. Can be 0, 1 or 2.
##' @param ... [internal] Only used by the generic method.
##' 
##' @return An object of class \code{calibrateType1}.
##' @seealso \code{link{autoplot.calibrateType1}} for a graphical display of the bias or of the type 1 error.
##' 
##' @author Brice Ozenne
##'
##' @examples
##' #### simulate data ####
##' m.Sim <- lvm(c(Y1[mu1:sigma]~1*eta,
##'                Y2[mu2:sigma]~1*eta,
##'                Y3[mu3:sigma]~1*eta,
##'                eta~beta1*Group+beta2*Gender))
##' latent(m.Sim) <- ~eta
##' categorical(m.Sim, labels = c("M","F")) <- ~Gender
##'
##' d <- lava::sim(m.Sim, 1e2)
##'
##' 
##' #### calibrate type 1 error on the estimated model ####
##' m <- lvm(Y1~eta,
##'          Y2~eta,
##'          Y3~eta,
##'          eta~Group+Gender)
##' e <- lava::estimate(m, data = d)
##' \dontshow{
##' res <- calibrateType1(e, null = "eta~Group", n.rep = 10)
##' }
##' \dontrun{
##' res <- calibrateType1(e, null = "eta~Group", n.rep = 100)
##' res <- calibrateType1(e, null = "eta~Group", n.rep = 100, cpus = 4)
##' }
##' summary(res)
##' 
##' @export
`calibrateType1` <-
  function(object, ...) UseMethod("calibrateType1")

## * calibrateType1.lvm
##' @rdname calibrateType1
##' @export
calibrateType1.lvm <- function(object, null, n, n.rep, F.test = FALSE, cluster = NULL,
                               generative.object = NULL, generative.coef = NULL, 
                               true.coef = NULL, n.true = 1e6, round.true = 2,              
                               bootstrap = FALSE, n.bootstrap = 1e3,
                               checkType1 = FALSE, checkType2 = FALSE,
                               dir.save = NULL, label.file = NULL,             
                               seed = NULL, cpus = 1, trace = 2, ...){

### ** test
    if(cpus>1){
        test.package <- try(requireNamespace("foreach"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'foreach\' \n",
                 "This package is necessary when argument \'cpus\' is greater than 1 \n")
        }

        if(cpus > parallel::detectCores()){
            stop("Argument \'cpus\' is greater than the number of available CPU cores \n",
                 "available CPU cores: ",parallel::detectCores(),"\n")
        }
    }
    
### ** prepare
    n.n <- length(n)

    ## *** generative model
    if(is.null(generative.object)){
        generative.object <- object
    }
    if(any(lava::manifest(object) %in% lava::manifest(generative.object) == FALSE)){
        missingVar <- lava::manifest(object)[lava::manifest(object) %in% lava::manifest(generative.object) == FALSE]
        stop("The object contains manifest variables that are not in the generative model \n",
             "missing manifest variables: \"",paste0(missingVar, collapse = "\" \""),"\" \n")
    }
    df.type_generative <- coefType(generative.object, as.lava = FALSE)
    name.param_generative <- df.type_generative[!is.na(df.type_generative$lava),"param"]
    n.param_generative <- length(name.param_generative)

    ## *** coef generative
    if(!is.null(generative.coef) && any(names(generative.coef) %in% name.param_generative == FALSE)){
        extraParam <- names(generative.coef)[names(generative.coef) %in% name.param_generative == FALSE]
        stop("Invalid argument \'generative.coef\': some of the coefficient names do not match those of the generative model \n",
             "extra coefficients: \"",paste0(extraParam, collapse = "\" \""),"\"\n")
    }

    ## *** coef of the fitted model
    if(is.null(true.coef)){
        if(trace>1){
            cat("  Estimate true coefficients using a sample size of n=",n.true," ", sep="")
        }
        e.true <- lava::estimate(object,
                                 data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE))
        coef.true <- coef(e.true)
        if(!is.null(round.true)){
            coef.true <- round(coef.true, digits = round.true)
        }
        if(trace>1){
            cat("- done \n")
        }
    }else{
        if(trace>1){
            cat("  Check true coefficients ")
        }
        n.true <- n[n.n]
        e.true <- lava::estimate(object, cluster = cluster,
                                 data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE))
        name.test <- names(coef(e.true))
        
        if(!identical(sort(name.test),sort(names(true.coef)))){
            extraNames <- setdiff(names(true.coef),name.test)
            missingNames <- setdiff(name.test,names(true.coef))
            stop("Names of the coefficients in argument \'true.coef\' does not matches those of the estimated coefficients \n",
                 "missing names: \"",paste0(missingNames, collapse = "\" \""),"\" \n",
                 "extra names: \"",paste0(extraNames, collapse = "\" \""),"\" \n")
        }
        
        coef.true <- true.coef
        if(trace>1){
            cat("- done \n")
        }
        
    }
    name.coef <- names(coef.true)
    n.coef <- length(name.coef)
       
    ## *** type of the coef of the fitted model
    df.type <- coefType(e.true, as.lava = FALSE)
    df.type <- df.type[df.type$name %in% name.coef,]
    type.coef <- setNames(df.type$detail, df.type$name)

    ## *** null hypothesis
    n.null <- length(null)
    if(any(null %in% name.coef == FALSE)){
        incorrect.name <- null[null %in% name.coef == FALSE]
        possible.name <- setdiff(name.coef, null)
        ls.name <- lapply(incorrect.name, function(iN){
            dist.tempo <- utils::adist(x = iN, y = possible.name)
            return(possible.name[which.min(dist.tempo)])
        })
        ex.name <- unique(unlist(ls.name))

        stop("Invalid argument \'null\': some of the coefficient names does not match those of the estimate model \n",
             "incorrect names: \"",paste(incorrect.name, collapse = "\" \""),"\" \n",
             "example of valid names: \"",paste(ex.name, collapse = "\" \""),"\"\n")
    }

    if(checkType1 && any(coef.true[null]!=0)){
        txtCoef <- paste(null[coef.true[null]!=0], collapse = "\" \"")
        stop("Control type 1 error: coefficients \"",txtCoef,"\" are not 0 while their belong to the null hypothesis\n")
    }
    if(checkType2 && any(coef.true[null]==0)){
        txtCoef <- paste(null[coef.true[null]==0], collapse = "\" \"")
        stop("Control type 2 error: coefficients \"",txtCoef,"\" are 0 while their belong to the null hypothesis\n")
    }

    res.C <- createContrast(null, name.param = name.coef, add.rowname = TRUE, rowname.rhs = FALSE)
    contrast <- res.C$contrast
    rhs <- res.C$null
    
### ** display
    if(trace>1){
        cat("  Settings: \n")
        cat("  > simulation for n=",paste(n,collapse = " "),"\n",sep="")
        cat("  > model: \n")
        print(object)
        cat("  > expected coefficients: \n")
        print(coef.true)
        if(bootstrap){
            cat("  > bootstrap: ",bootstrap,"\n")
        }
        if(!is.null(seed)){
            cat("  > seed: ",seed,"\n")
        }        
    }

    
### ** loop
    store.coef <- null
    if(F.test){
        store.coef <- c(store.coef, "global")
    }
    n.store <- length(store.coef)
    
    if(trace>1){
        cat("\n")
        cat(" Perform simulation: \n")
    }

    if(cpus>1){

        ## *** define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.rep, style = 3)          
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        ## *** link to foreach
        doParallel::registerDoParallel(cl)

        ## *** export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(requireNamespace("lava", quietly = TRUE))
            suppressPackageStartupMessages(requireNamespace("lavaSearch2", quietly = TRUE))
        })
        
        ## *** seed
        cpus.name <- unlist(parallel::clusterApply(cl = cl, 1:cpus, function(x){
            myName <- paste(Sys.info()[["nodename"]], Sys.getpid(), sep="-")
            return(myName)
        }))
        if(length(seed)==0){
            seed <- rep(NA, cpus)
        }else if(length(seed)==1){
            set.seed(seed)
            seed <- sample(1:1e5,size=cpus,replace=FALSE)                
        }else{
            if(length(seed)!=cpus){
                stop("Length of argument \'seed\' does not match argument \'cpus\' \n")
            }            
        }
        names(seed) <- cpus.name

        
        ## *** parallel computation
        toExport <- c(".warperType1", "cpus.name")

        iRep <- NULL # [:for CRAN check] foreach
        resSim <- foreach::`%dopar%`(
                               foreach::foreach(iRep = 1:n.rep,
                                                .export = toExport),{ # iRep <- 1

                                                    if(trace>0){utils::setTxtProgressBar(pb, iRep)}

                                                    myName <- paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
                                                    iSeed <- seed[myName]
                                                    
                                                    ls.pvalue <- vector(mode = "list", length = n.n)
                                                    ls.bias <- vector(mode = "list", length = n.n)
                                                    iIndex <- 1
                                                    
                                                    for(iN in n){ # iN <- n[1]
                                                        iRes <- .warperType1(iRep,
                                                                             n = iN,
                                                                             generative.object = generative.object,
                                                                             generative.coef = generative.coef,
                                                                             object = object,
                                                                             cluster = cluster,
                                                                             coef.true = coef.true,
                                                                             type.coef = type.coef,
                                                                             name.coef = name.coef,
                                                                             store.coef = store.coef,
                                                                             n.coef = n.coef,
                                                                             n.store = n.store,
                                                                             F.test = F.test,
                                                                             null = null,
                                                                             contrast = contrast,
                                                                             rhs = rhs,
                                                                             bootstrap = bootstrap,
                                                                             n.bootstrap = n.bootstrap,
                                                                             seed = iSeed)

                                                        ls.pvalue[[iIndex]] <- iRes$pvalue
                                                        ls.bias[[iIndex]] <- iRes$bias
                                                        iIndex <- iIndex + 1
                                                    }
                                                    
                                                    return(list(pvalue = do.call("rbind",ls.pvalue),
                                                                bias = do.call("rbind",ls.bias)))
                                                })
    
        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
        
        ## *** post process
        dt.pvalue <- do.call("rbind",lapply(resSim,"[[","pvalue"))
        dt.pvalue <- dt.pvalue[order(dt.pvalue$n,dt.pvalue$rep),,drop=FALSE]
        dt.bias <- do.call("rbind",lapply(resSim,"[[","bias"))
        dt.bias <- dt.bias[order(dt.bias$n,dt.bias$rep),,drop=FALSE]
        
    }else{

        ## *** filename
        if(is.null(label.file)){label.file <- seed}
        filename_tempo.pvalue <- paste0("type1error-S",label.file,"(tempo).rds")
        filename_tempo.bias <- paste0("bias-S",label.file,"(tempo).rds")
        filename.pvalue <- gsub("\\(tempo\\)","",filename_tempo.pvalue)
        filename.bias <- gsub("\\(tempo\\)","",filename_tempo.bias)

        if(!is.null(dir.save)){
            validPath(dir.save, type = "dir")
        }

        if(!is.null(dir.save)){
            cat("  > export results in ",dir.save,"\n")
        }

        if(trace>0){
            test.package <- try(requireNamespace("pbapply"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'pbapply\' \n",
                     "This package is necessary when argument \'trace\' is TRUE \n")
            }
            FCTapply <- pbapply::pblapply
        }else{
            FCTapply <- lapply
        }
        if(!is.null(seed)){
            set.seed(seed)
        }else{
            seed <- NA
        }
        
        ## *** sequential simulation
        dt.pvalue <- NULL
        dt.bias <- NULL

        for(iN in n){ ## iN <- n[1]

            if(trace>0){
                cat("  > sample size=",iN,"\n", sep = "")                    
            }

            resSim <- do.call(FCTapply, args = list(X = 1:n.rep, FUN = function(iRep){
                .warperType1(iRep,
                             n = iN,
                             generative.object = generative.object,
                             generative.coef = generative.coef,
                             object = object, cluster = cluster,                             
                             coef.true = coef.true, type.coef = type.coef, name.coef = name.coef,
                             store.coef = store.coef, n.coef = n.coef, n.store = n.store,
                             F.test = F.test, null = null, contrast = contrast, rhs = rhs,
                             bootstrap = bootstrap,
                             n.bootstrap = n.bootstrap,
                             seed = seed)
            }))
            dt.pvalue <- rbind(dt.pvalue,
                               do.call("rbind",lapply(resSim,"[[","pvalue"))
                               )
            dt.bias <- rbind(dt.bias,
                             do.call("rbind",lapply(resSim,"[[","bias")))
                             
            
            ## export (tempo)
            if(!is.null(dir.save)){
                saveRDS(dt.pvalue, file = file.path(dir.save,filename_tempo.pvalue))
                saveRDS(dt.bias, file = file.path(dir.save,filename_tempo.bias))
            }
            
        }   
    }

    ## ** export
    if(!is.null(dir.save)){
        saveRDS(dt.pvalue, file = file.path(dir.save,filename.pvalue))
        saveRDS(dt.bias, file = file.path(dir.save,filename.bias))
    }
    out <- list(p.value = dt.pvalue,
                bias = dt.bias,
                e.true = e.true,
                null = null)
    class(out) <- append("calibrateType1",class(out))
    return(out)


}

## * calibrateType1.lvmfit
##' @rdname calibrateType1
##' @export
calibrateType1.lvmfit <- function(object, null, n.rep, F.test = FALSE,
                                  bootstrap = FALSE, n.bootstrap = 1e3,
                                  seed = NULL, trace = 2, cpus = 1, ...){

    ## ** Prepare
    ## *** model
    object.model <- object$model

    ## *** coef
    coef.true <- coef(object)
    name.coef <- names(coef.true)
    if(any(null %in% name.coef == FALSE)){
        txt <- null[null %in% name.coef == FALSE]
        txt2 <- setdiff(name.coef, null)
        stop("Argument \'null\' does not match the names of the model coefficients \n",
             "Incorrect null: \"",paste(txt, collapse = "\" \""),"\" \n",
             "Possible null: \"",paste(txt2, collapse = "\" \""),"\" \n")
    }
    coef.true[null] <- 0

    ## *** data
    n <- object$data$n

    ## ** Run
    out <- calibrateType1(object.model,
                          null = null,
                          n = n,
                          n.rep = n.rep,
                          F.test = F.test,
                          generative.object = object.model,
                          generative.coef = coef.true, 
                          true.coef = coef.true,              
                          bootstrap = bootstrap,
                          n.bootstrap = n.bootstrap,
                          checkType1 = FALSE,
                          checkType2 = FALSE,
                          dir.save = NULL,
                          label.file = NULL,             
                          seed = seed,
                          cpus = cpus,
                          trace = trace)


    ## ** Export
    return(out)
}

## * .warperType1
.warperType1 <- function(iRep, n, generative.object, generative.coef,
                         object, cluster,
                         coef.true, type.coef, name.coef, store.coef, n.coef, n.store,
                         F.test, null, contrast, rhs,
                         bootstrap, n.bootstrap,
                         seed){

    ls.iP <- list()  ## temporary
    out <- list()
    
    ## ** simulation
    dt.sim <- lava::sim(generative.object, n = n, p = generative.coef, latent = FALSE)

    ## ** model adjustement
    e.lvm <- lava::estimate(object, data = dt.sim, cluster = cluster)
    if(e.lvm$opt$convergence==1){return(list(pvalue=NULL,bias=NULL))} ## exclude lvm that has not converged
    if(any(eigen(getVarCov2(e.lvm))$values<=0)){return(list(pvalue=NULL,bias=NULL))} ## exclude lvm where the residual covariance matrix is not semipositive definite

    e.lvm.Satt <- e.lvm
    testError.Satt <- try(sCorrect(e.lvm.Satt) <- FALSE, silent = TRUE)
    e.lvm.KR <- e.lvm
    testError.KR <- try(suppressWarnings(sCorrect(e.lvm.KR, safeMode = TRUE) <- TRUE), silent = TRUE)
            
    ## ** coefficients
    coef.original <- coef(e.lvm)
    if(!inherits(testError.KR,"try-error")){
        ## check whether adjusted residuals could be computed (otherwise adjust.n=FALSE)
        test.warning <- inherits(attr(e.lvm.KR$sCorrect,"warning"),"try-error")
        niter.correct <- e.lvm.KR$sCorrect$opt$iterations
        coef.corrected <- e.lvm.KR$sCorrect$param
    }else{
        niter.correct <- NA
        test.warning <- NA
        coef.corrected <- setNames(rep(NA,n.coef),name.coef)
    }
            
    ## ** no correction
    ## get Wald tests
    eS.ML <- try(summary(e.lvm)$coef[,c("Estimate","P-value")], silent = TRUE)
    if("try-error" %in% class(eS.ML)){return(list(pvalue=NULL,bias=NULL))} ## exclude lvm where we cannot compute the summary
    if(F.test){
        F.ML <- lava::compare(e.lvm, par = null)
        eS.ML <- rbind(eS.ML, global = c(Estimate = F.ML$statistic, "P-value" = F.ML$p.value))
    }

    if(is.null(cluster)){
        if(!inherits(testError.Satt,"try-error")){                
            eS.robustML <- compare2(e.lvm.Satt, robust = TRUE, df = FALSE,
                                    contrast = contrast, null = rhs,
                                    F.test = F.test, as.lava = FALSE)[,c("estimate","p-value")]
            names(eS.robustML) <- c("Estimate","P-value")
        }else{
            eS.robustML <- try(estimate(e.lvm)$coefmat[,c("Estimate","P-value")], silent = TRUE)
            if(inherits(eS.robustML,"try-error")){
                eS.robustML <- matrix(as.numeric(NA), ncol = 2, nrow = n.coef, dimnames = list(name.coef, c("Estimate","P-value")))
            }

            if(F.test){
                eS.robustML <- rbind(eS.robustML, global = rep(NA,2))
            }
        }
    }else{
        eS.robustML <- eS.ML
        eS.ML[] <- NA
    }

    ## store results
    ls.iP$p.Ztest <- eS.ML[store.coef,"P-value"]
    ls.iP$p.robustZtest <-  eS.robustML[store.coef,"P-value"]
            
    ## ** Sattterwaith correction
    if(!inherits(testError.Satt,"try-error")){
        ## get Wald tests
        eS.Satt <- compare2(e.lvm.Satt, robust = FALSE, df = TRUE,
                            contrast = contrast, null = rhs,
                            F.test = F.test, as.lava = FALSE)
        eS.robustSatt <- compare2(e.lvm.Satt, robust = TRUE, df = TRUE,
                                  contrast = contrast, null = rhs,
                                  F.test = F.test, as.lava = FALSE)

        ## store results
        ls.iP$p.Satt <- eS.Satt[store.coef,"p-value"]
        ls.iP$p.robustSatt <- eS.robustSatt[store.coef,"p-value"]
    }else{
        ls.iP$p.Satt <- rep(as.numeric(NA), n.store)
        ls.iP$p.robustSatt <- rep(as.numeric(NA), n.store)
    }

    ## ** small sample correction
    if(!inherits(testError.KR,"try-error")){
        ## get Wald tests
        eS.SSC <- compare2(e.lvm.KR, robust = FALSE, df = FALSE,
                           contrast = contrast, null = rhs,
                           F.test = F.test, as.lava = FALSE)
        eS.robustSSC <- compare2(e.lvm.KR, robust = TRUE, df = FALSE,
                                 contrast = contrast, null = rhs,
                                 F.test = F.test, as.lava = FALSE)
        ## store results
        ls.iP$p.SSC <- eS.SSC[store.coef,"p-value"]
        ls.iP$p.robustSSC <- eS.robustSSC[store.coef,"p-value"]
    }else{
        ls.iP$p.SSC <- rep(as.numeric(NA), n.store)
        ls.iP$p.robustSSC <- rep(as.numeric(NA), n.store)
    }

    ## ** Sattterwaith correction with small sample correction
    if(!inherits(testError.KR,"try-error")){
        ## get Wald tests
        eS.KR <- compare2(e.lvm.KR, robust = FALSE, df = TRUE,
                          contrast = contrast, null = rhs,
                          F.test = F.test, as.lava = FALSE)
        eS.robustKR <- compare2(e.lvm.KR, robust = TRUE, df = TRUE,
                                contrast = contrast, null = rhs,
                                F.test = F.test, as.lava = FALSE)
                
        ## store results
        ls.iP$p.KR <- eS.KR[store.coef,"p-value"]
        ls.iP$p.robustKR <- eS.robustKR[store.coef,"p-value"]
    }else{
        ls.iP$p.KR <- rep(as.numeric(NA), n.store)
        ls.iP$p.robustKR <- rep(as.numeric(NA), n.store)
    }
    ## ** bootstrap
    if(bootstrap>0){
        e.boot <- eval(parse(text = "butils::bootReg(e.lvm, type = \"coef\", n.boot = n.bootstrap"))

        index.coef.boot <- match(null, name.coef)
        boot.perc <- summary(e.boot, p.value = TRUE, type = "perc", print = FALSE, index = index.coef.boot)
        boot.stud <- summary(e.boot, p.value = TRUE, type = "stud", print = FALSE, index = index.coef.boot)
        boot.bca <- summary(e.boot, p.value = TRUE, type = "bca", print = FALSE, index = index.coef.boot)

        ls.iP$p.bootPerc <- boot.perc[store.coef,"p.value"]
        ls.iP$p.bootStud <- boot.stud[store.coef,"p.value"]
        ls.iP$p.bootBca <- boot.bca[store.coef,"p.value"]
    }

    ## ** metainformation
    out$pvalue <- cbind(data.frame(n = n,
                                   rep = iRep,
                                   seed = seed,
                                   nboot = n.bootstrap,
                                   niter = niter.correct,
                                   warning = test.warning,
                                   link = store.coef,
                                   stringsAsFactors = FALSE),
                        do.call(cbind,ls.iP))
    rownames(out$pvalue) <- NULL
    
    ## ** collect result (bias)
    out$bias <- data.frame(n = n,
                           rep = iRep,
                           seed = seed,
                           niter = niter.correct,
                           warning = test.warning,
                           estimate.truth = as.double(coef.true),
                           estimate.ML = as.double(coef.original[name.coef]),
                           estimate.MLcorrected = as.double(coef.corrected[name.coef]),
                           name = names(coef.true),
                           type = type.coef[name.coef],
                           stringsAsFactors = FALSE)
    rownames(out$bias) <- NULL

    ## ** export
    return(out)
}
    
######################################################################
### calibrateType1.R ends here

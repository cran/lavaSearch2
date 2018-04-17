### calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (10:23) 
## Version: 
## Last-Updated: apr 17 2018 (10:55) 
##           By: Brice Ozenne
##     Update #: 332
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
##' @name calibrateType
##' 
##' @param object a \code{lvm} object defining the model to be fitted.
##' @param null [character vector] names of the coefficient whose value will be tested against 0. 
##' @param n [integer vector, >0] sample size(s) considered in the simulation study.
##' @param n.rep [integer, >0] number of simulations per sample size.
##' @param generative.object [lvm] object defining the statistical model generating the data.
##' @param generative.coef [name numeric vector] values for the parameters of the generative model.
##' Can also be \code{NULL}: in such a case the coefficients are set to default values decided by lava (usually 0 or 1).
##' @param true.coef [name numeric vector] expected values for the parameters of the fitted model.
##' @param n.true [integer, >0] sample size at which the estimated coefficients will be a reliable approximation of the true coefficients.
##' @param round.true [integer, >0] the number of decimal places to be used for the true value of the coefficients. No rounding is done if \code{NULL}.
##' @param bootstrap [logical] should bootstrap resampling be performed?
##' @param type.bootstrap [character vector]
##' @param n.bootstrap [integer, >0] the number of bootstrap sample to be used for each bootstrap.
##' @param checkType1 [logical] returns an error if the coefficients associated to the null hypotheses do not equal 0.
##' @param checkType2 [logical] returns an error if the coefficients associated to the null hypotheses equal 0.
##' @param dir.save [character] path to the directory were the results should be exported.
##' Can also be \code{NULL}: in such a case the results are not exported.
##' @param F.test [logical] should a multivariate Wald test be perform testing simultaneously all the null hypotheses?
##' @param label.file [character] element to include in the file name.
##' @param seed [integer, >0] seed value that will be set at the beginning of the simulation to enable eproducibility of the results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param trace [interger] should the execution of the function be trace. Can be 0, 1 or 2.
##' 
##' @return An object of class \code{calibrateType1}.
##' @seealso \code{link{autoplot.calibrateType1}} for a graphical display of the bias or of the type 1 error.
##' 
##' @author Brice Ozenne
##'
##' @examples
##' #### generative model ####
##' m.Sim <- lvm(c(Y1[mu1:sigma]~1*eta,
##'                Y2[mu2:sigma]~1*eta,
##'                Y3[mu3:sigma]~1*eta,
##'                eta~beta1*Group+beta2*Gender))
##' latent(m.Sim) <- ~eta
##' categorical(m.Sim, labels = c("M","F")) <- ~Gender
##'
##' vec.par <- c(mu2 = 0, ## mu1 is set to 0 by default
##'              mu3 = -0.37,
##'              "eta" = -2.74,
##'              beta1 = -1.20,
##'              beta2 = 0,
##'              sigma = 1.46,
##'              "eta~~eta" = 1.63)
##'
##' #### parameters to test ####
##' null <- c("Y2","eta~GenderF")
##'
##' #### launch simulation: same model generation fit ####
##' \dontrun{
##' res <- calibrateType1(m.Sim, null = null, n = c(20,30,40), n.rep = 50, generative.coef = vec.par)
##' autoplot(res, type = "bias")
##' autoplot(res, type = "type1error")
##' }
##' \dontshow{
##' res <- calibrateType1(m.Sim, null = null, n = c(20,30,40), n.true = 1e3, n.rep = 2, generative.coef = vec.par)
##' }

## * calibrateType1
##' @rdname calibrateType1
##' @export
calibrateType1 <- function(object, null, n, n.rep, F.test = FALSE,
                           generative.object = NULL, generative.coef = NULL, 
                           true.coef = NULL, n.true = 1e6, round.true = 2,              
                           bootstrap = FALSE, type.bootstrap = c("perc","stud","bca"), n.bootstrap = 1e3,
                           checkType1 = FALSE, checkType2 = FALSE,
                           dir.save = NULL, label.file = NULL,             
                           seed = NULL, trace = 2){

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
            cat("* estimate true coefficients using a sample size of n=",n.true," ", sep="")
        }
        e.true <- lava::estimate(object, data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE))
        coef.true <- coef(e.true)
        if(!is.null(round.true)){
            coef.true <- round(coef.true, digits = round.true)
        }
        if(trace>1){
            cat("- done \n")
        }
    }else{
        if(trace>1){
            cat("* check true coefficients ")
        }
        n.true <- n[n.n]
        e.true <- lava::estimate(object, data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE))
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
        stop("Invalid argument \'null\': some of the coefficient names does not match those of the estimate model \n",
             "incorrect names: \"",paste(incorrect.name, collapse = "\" \""),"\" \n")
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
    
    ## *** filename
    if(is.null(label.file)){label.file <- seed}
    filename_tempo.pvalue <- paste0("type1error-S",label.file,"(tempo).rds")
    filename_tempo.bias <- paste0("bias-S",label.file,"(tempo).rds")
    filename.pvalue <- gsub("\\(tempo\\)","",filename_tempo.pvalue)
    filename.bias <- gsub("\\(tempo\\)","",filename_tempo.bias)

    if(!is.null(dir.save)){
        validPath(dir.save, type = "dir")
    }

### ** display
    if(trace>1){
        cat("* settings: \n")
        cat("  > simulation for n=",paste(n,collapse = " "),"\n",sep="")
        cat("  > model: \n")
        print(object)
        cat("  > expected coefficients: \n")
        print(coef.true)
        cat("  > bootstrap: ",bootstrap,"\n")
        if(!is.null(seed)){
            cat("  > seed: ",seed,"\n")
        }
        if(!is.null(dir.save)){
            cat("  > export results in ",dir.save,"\n")
        }
        
    }
    
### ** loop
    dt.pvalue <- NULL
    dt.bias <- NULL
    store.coef <- null
    if(F.test){
        store.coef <- c(store.coef, "global")
    }
    n.store <- length(store.coef)
    if(!is.null(seed)){
        set.seed(seed)
    }else{
        seed <- NA
    }

    if(trace>1){cat("* perform simulation: \n")}
    for(iN in 1:n.n){

        if(trace>0){cat("  > sample size=",n[iN],": ", sep = "")}
        n.tempo <- n[iN]

        for(iRep in 1:n.rep){
            if(trace>0){cat(iRep," ")}
            ls.iP <- list()
            
            ## *** simulation
            dt.sim <- lava::sim(generative.object, n = n.tempo, p = generative.coef, latent = FALSE)

            ## *** model adjustement
            e.lvm <- lava::estimate(object, data = dt.sim)
            if(e.lvm$opt$convergence==1){next} ## exclude lvm that has not converged
            if(any(eigen(getVarCov2(e.lvm))$values<=0)){next} ## exclude lvm where the residual covariance matrix is not semipositive definite

            e.lvm.Satt <- e.lvm
            testError.Satt <- try(sCorrect(e.lvm.Satt) <- FALSE)
            e.lvm.KR <- e.lvm
            testError.KR <- try(suppressWarnings(sCorrect(e.lvm.KR, safeMode = TRUE) <- TRUE))
            
            ## *** coefficients
            coef.original <- coef(e.lvm)
            if(!inherits(testError.KR,"try-error")){
                ## check whether adjusted residuals could be computed (otherwise adjust.n=FALSE)
                test.warning <- inherits(attr(e.lvm.KR$sCorrect,"warning"),"try-error")

                coef.corrected <- e.lvm.KR$sCorrect$param
            }else{
                coef.corrected <- setNames(rep(NA,n.coef),name.coef)
                test.warning <- NA
            }
            
            ## *** no correction
            ## get Wald tests
            eS.ML <- try(summary(e.lvm)$coef[,c("Estimate","P-value")],silent = TRUE)
            if("try-error" %in% class(eS.ML)){next} ## exclude lvm where we cannot compute the summary
            if(F.test){
                F.ML <- lava::compare(e.lvm, par = null)
                eS.ML <- rbind(eS.ML, global = c(Estimate = F.ML$statistic, "P-value" = F.ML$p.value))
            }
            
            if(!inherits(testError.Satt,"try-error")){                
                eS.robustML <- compare2(e.lvm.Satt, robust = TRUE, df = FALSE,
                                        contrast = contrast, null = rhs,
                                        F.test = F.test, as.lava = FALSE)[,c("estimate","p-value")]
                names(eS.robustML) <- c("Estimate","P-value")
            }else{
                eS.robustML <- try(estimate(e.lvm)$coefmat[,c("Estimate","P-value")],silent = TRUE)
                if(inherits(eS.robustML,"try-error")){
                    eS.robustML <- matrix(as.numeric(NA), ncol = 2, nrow = n.coef, dimnames = list(name.coef, c("Estimate","P-value")))
                }
                if(F.test){
                    eS.robustML <- rbind(eS.robustML, global = rep(NA,2))
                }
            }

            ## store results
            ls.iP$p.Ztest <- eS.ML[store.coef,"P-value"]
            ls.iP$p.robustZtest <-  eS.robustML[store.coef,"P-value"]
            
            ## *** Sattterwaith correction
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

            ## *** small sample correction
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

            ## *** Sattterwaith correction with small sample correction
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
            ## *** bootstrap
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

            ## *** metainformation
            iDT.pvalue <- cbind(data.frame(n = n.tempo,
                                           rep = iRep,
                                           seed = seed,
                                           nboot = n.bootstrap,
                                           niter = e.lvm.KR$sCorrect$opt$iterations,
                                           warning = test.warning,
                                           link = store.coef,
                                           stringsAsFactors = FALSE),
                                do.call(cbind,ls.iP))
            rownames(iDT.pvalue) <- NULL
            dt.pvalue <- rbind(dt.pvalue, iDT.pvalue)

            ## *** collect result (bias)
            iDT.bias <- data.frame(n = n.tempo,
                                   rep = iRep,
                                   seed = seed,
                                   niter = e.lvm.KR$sCorrect$opt$iterations,
                                   warning = test.warning,
                                   estimate.truth = as.double(coef.true),
                                   estimate.ML = as.double(coef.original[name.coef]),
                                   estimate.MLcorrected = as.double(coef.corrected[name.coef]),
                                   name = names(coef.true),
                                   type = type.coef[name.coef],
                                   stringsAsFactors = FALSE)
            rownames(iDT.bias) <- NULL
            dt.bias <- rbind(dt.bias, iDT.bias)
        }

        ## *** export (tempo)
        if(!is.null(dir.save)){
            saveRDS(dt.pvalue, file = file.path(dir.save,filename_tempo.pvalue))
            saveRDS(dt.bias, file = file.path(dir.save,filename_tempo.bias))
        }
        if(trace>0){cat("\n")}
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


######################################################################
### calibrateType1.R ends here

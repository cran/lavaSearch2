### iidJack.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (09:15) 
## Version: 
## last-updated: feb  5 2018 (18:15) 
##           By: Brice Ozenne
##     Update #: 303
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - iidJack
#' @title Jackknife iid Decomposition from Model Object
#' @description Extract iid decomposition (i.e. influence function) from model object.
#'
#' @name iidJack
#' 
#' @param object a object containing the model.
#' @param data [data.frame] dataset used to perform the jackknife.
#' @param grouping [vector] variable defining cluster of observations that will be simultaneously removed by the jackknife.
#' @param ncpus [integer >0] the number of processors to use.
#' If greater than 1, the fit of the model and the computation of the influence function for each jackknife sample is performed in parallel. 
#' @param keep.warnings [logical] keep warning messages obtained when estimating the model with the jackknife samples.
#' @param keep.error [logical]keep error messages obtained when estimating the model with the jackknife samples.
#' @param init.cpus [logical] should the processors for the parallel computation be initialized?
#' @param trace [logical] should a progress bar be used to trace the execution of the function
#' @param ... [internal] only used by the generic method.
#'
#' @return A matrix with in row the samples and in columns the parameters.
#' 
#' @examples
#' n <- 20
#'
#' #### glm ####
#' set.seed(10)
#' m <- lvm(y~x+z)
#' distribution(m, ~y+z) <- binomial.lvm("logit")
#' d <- lava::sim(m,n)
#' g <- glm(y~x+z,data=d,family="binomial")
#' iid1 <- iidJack(g, ncpus = 1)
#' iid2 <- lava::iid(g)
#' quantile(iid1-iid2)
#' vcov(g)
#' colSums(iid2^2)
#' colSums(iid1^2)
#' 
#' #### Cox model ####
#' library(survival)
#' data(Melanoma, package = "riskRegression")
#' m <- coxph(Surv(time,status==1)~ici+age, data = Melanoma, x = TRUE, y = TRUE)
#' 
#' \dontrun{
#' ## require riskRegression > 1.4.3
#' if(utils::packageVersion("riskRegression") > "1.4.3"){
#' library(riskRegression)
#' iid1 <- iidJack(m)
#' iid2 <- iidCox(m)$IFbeta
#'
#' apply(iid1,2,sd)
#'
#' print(iid2)
#' 
#' apply(iid2,2,sd)
#'   }
#' }
#' 
#' #### LVM ####
#' set.seed(10)
#'
#' mSim <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ 1*eta)
#' latent(mSim) <- ~eta
#' categorical(mSim, K=2) <- ~G
#' transform(mSim, Id ~ eta) <- function(x){1:NROW(x)}
#' dW <- lava::sim(mSim, n, latent = FALSE)
#' dL <- reshape2::melt(dW, id.vars = c("G","Id"),
#'                      variable.name = "time", value.name = "Y")
#' dL$time <- gsub("Y","",dL$time)
#'
#' m1 <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ 1*eta)
#' latent(m1) <- ~eta
#' regression(m1) <- eta ~ G
#' e <- estimate(m1, data = dW)
#' \dontrun{
#' iid1 <- iidJack(e)
#' iid2 <- iid(e)
#' attr(iid2, "bread") <- NULL
#'
#' apply(iid1,2,sd)
#' apply(iid2,2,sd)
#' quantile(iid2 - iid1)
#' }
#' 
#' library(nlme)
#' e2 <- lme(Y~G+time, random = ~1|Id, weights = varIdent(form =~ 1|Id), data = dL)
#' e2 <- lme(Y~G, random = ~1|Id, data = dL)
#' \dontrun{
#' iid3 <- iidJack(e2)
#' apply(iid3,2,sd)
#' }
#'
#' @concept iid decomposition
#' @export
iidJack <- function(object,...) UseMethod("iidJack")

## * method iidJack.default
#' @rdname iidJack
#' @export
iidJack.default <- function(object,data=NULL,grouping=NULL,ncpus=1,
                            keep.warnings=TRUE, keep.error=TRUE,
                            init.cpus=TRUE,trace=TRUE,...) {
    
    estimate.lvm <- lava_estimate.lvm

    ## ** extract data
    if(is.null(data)){
        myData <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE)
    }else{ 
        myData <-  as.data.frame(data)
    }
    
    n.obs <- NROW(myData)
    if(any(class(object) %in% "lme")){
        getCoef <- nlme::fixef
    }else{
        getCoef <- coef
    }
    coef.x <- getCoef(object)
    names.coef <- names(coef.x)
    n.coef <- length(coef.x)

    ## ** update formula/model when defined by a variable and not in the current namespace
    if(length(object$call[[2]])==1){
        modelName <- as.character(object$call[[2]])
        if(modelName %in% ls() == FALSE){
            assign(modelName, value = evalInParentEnv(object$call[[2]]))
        }
    }

    ## ** define the grouping level for the data
    if(is.null(grouping)){
        if(any(class(object)%in%c("lme","gls","nlme"))){
            myData$XXXgroupingXXX <- as.vector(apply(object$groups,2,interaction))
        }else{
            myData$XXXgroupingXXX <- 1:n.obs
        }
        grouping <- "XXXgroupingXXX"        
    }else{
        if(length(grouping)>1){
            stop("grouping must refer to only one variable \n")
        }
        if(grouping %in% names(myData) == FALSE){
            stop("variable defined in grouping not found in data \n")
        }
    }
    myData[,grouping] <- as.character(myData[,grouping])
    Ugrouping <- unique(myData[,grouping])
    n.group <- length(Ugrouping)

    ## ** warper
    warper <- function(i){ # i <- "31"
        newData <- subset(myData, subset = myData[[grouping]]!=i)
        xnew <- tryWithWarnings(stats::update(object, data = newData))
        if(!is.null(xnew$error)){
            xnew$value <- rep(NA, n.coef)
        }else{
            xnew$value <- getCoef(xnew$value)
        }
        return(xnew)
    }
    # warper("31")
    
    ## ** parallel computations: get jackknife coef
    if(ncpus>1){
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

        estimator <- as.character(object$call[[1]]) 

        vec.packages <- c("lava")
        possiblePackage <- gsub("package:","",utils::getAnywhere(estimator)$where[1])
        existingPackage <- as.character(utils::installed.packages()[,"Package"])

        ls.call <- as.list(object$call)
        test.length <- which(unlist(lapply(ls.call, length))==1)
        test.class <- which(unlist(lapply(ls.call, function(cc){
            (class(c) %in% c("numeric","character","logical")) == FALSE
        })))
        test.class <- which(unlist(lapply(ls.call, class)) %in% c("numeric","character","logical") == FALSE)
    
        indexExport <- intersect(test.class,test.length)
        toExport <- sapply(ls.call[indexExport], as.character)
    
        if(possiblePackage %in% existingPackage){
            vec.packages <- c(vec.packages,possiblePackage)
        }
        if(length(object$call$data)==1){
            toExport <- c(toExport,as.character(object$call$data))
        }
        if(length(object$call$formula)==1){
            toExport <- c(toExport,as.character(object$call$formula))
        }
        if(length(object$call$fixed)==1){
            toExport <- c(toExport,as.character(object$call$fixed))        
        }
        toExport <- c(unique(toExport),"tryWithWarnings")

        #sapply(as.list(object$call),as.character)
        i <- NULL # [:for CRAN check] foreach
        resLoop <- foreach::`%dopar%`(
                                foreach::foreach(i = 1:n.group, .packages =  vec.packages,
                                                 .export = toExport),{                                                      
                                                     if(trace){
                                                         if(!exists("pb")){
                                                             pb <- tcltk::tkProgressBar("iidJack:", min=1, max=n.group)
                                                         }
                                                         tcltk::setTkProgressBar(pb, i)
                                                     }
                                                     warper(Ugrouping[i])
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
            resLoop <- pbapply::pblapply(Ugrouping, warper)
            
        }else{
            resLoop <- lapply(Ugrouping, warper)
        }
        

    }
    coefJack <- do.call(rbind, lapply(resLoop,"[[","value"))
    rownames(coefJack) <- 1:n.group

    ## ** post treatment: from jackknife coef to the iid decomposition
    # defined as (n-1)*(coef-coef(-i))
    # division by n to match output of lava, i.e. IF/n
    iidJack <- -(n.group-1)/n.group * sweep(coefJack, MARGIN = 2, STATS = coef.x, FUN = "-")
    colnames(iidJack) <- names.coef

    if(keep.warnings){
        ls.warnings <- lapply(resLoop,"[[","warnings")
        names(ls.warnings) <- 1:n.group
        ls.warnings <- Filter(Negate(is.null), ls.warnings)
        attr(iidJack,"warnings") <- ls.warnings
    }
    if(keep.error){
        ls.error <- lapply(resLoop,"[[","error")
        names(ls.error) <- 1:n.group
        ls.error <- Filter(Negate(is.null), ls.error)
        attr(iidJack,"error") <- ls.error
    }
    return(iidJack)
}
    
#----------------------------------------------------------------------
### iidJack.R ends here

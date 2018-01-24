### iidJack.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (09:15) 
## Version: 
## last-updated: jan 19 2018 (10:23) 
##           By: Brice Ozenne
##     Update #: 283
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - iidJack
#' @title Jacknife iid Decomposition from Model Object
#' @description Extract iid decomposition (i.e. influence function) from model object.
#'
#' @name iidJack
#' 
#' @param x model object.
#' @param data dataset used to perform the jacknife.
#' @param grouping variable defining cluster of observations that will be simultaneously removed by the jackknife.
#' @param ncpus number of cpus available for parallel computation.
#' @param keep.warnings keep warning messages obtained when estimating the model with the jackknife samples.
#' @param keep.error keep error messages obtained when estimating the model with the jackknife samples.
#' @param initCpus should the parallel computation be initialized?
#' @param trace should a progress bar be used to trace the execution of the function
#' @param ... additional arguments.
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
#' @export
iidJack <- function(x,...) UseMethod("iidJack")

## * method iidJack.default
#' @rdname iidJack
#' @export
iidJack.default <- function(x,data=NULL,grouping=NULL,ncpus=1,
                            keep.warnings=TRUE, keep.error=TRUE,
                            initCpus=TRUE,trace=TRUE,...) {
    
    estimate.lvm <- lava_estimate.lvm

    ## ** extract data
    if(is.null(data)){
        myData <- extractData(x, design.matrix = FALSE, as.data.frame = TRUE)
    }else{ 
        myData <-  as.data.frame(data)
    }
    
    n.obs <- NROW(myData)
    if(any(class(x) %in% "lme")){
        getCoef <- nlme::fixef
    }else{
        getCoef <- coef
    }
    coef.x <- getCoef(x)
    names.coef <- names(coef.x)
    n.coef <- length(coef.x)

    ## ** update formula/model when defined by a variable and not in the current namespace
    if(length(x$call[[2]])==1){
        modelName <- as.character(x$call[[2]])
        if(modelName %in% ls() == FALSE){
            assign(modelName, value = evalInParentEnv(x$call[[2]]))
        }
    }

    ## ** define the grouping level for the data
    if(is.null(grouping)){
        if(any(class(x)%in%c("lme","gls","nlme"))){
            myData$XXXgroupingXXX <- as.vector(apply(x$groups,2,interaction))
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
        xnew <- tryWithWarnings(stats::update(x, data = newData))
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
        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doParallel::registerDoParallel(cl)
        }
 
        if(trace > 0){
            parallel::clusterExport(cl, varlist = "trace")
        }

        estimator <- as.character(x$call[[1]]) 

        vec.packages <- c("lava")
        possiblePackage <- gsub("package:","",utils::getAnywhere(estimator)$where[1])
        existingPackage <- as.character(utils::installed.packages()[,"Package"])

        ls.call <- as.list(x$call)
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
        if(length(x$call$data)==1){
            toExport <- c(toExport,as.character(x$call$data))
        }
        if(length(x$call$formula)==1){
            toExport <- c(toExport,as.character(x$call$formula))
        }
        if(length(x$call$fixed)==1){
            toExport <- c(toExport,as.character(x$call$fixed))        
        }
        toExport <- c(unique(toExport),"tryWithWarnings")

        #sapply(as.list(x$call),as.character)
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
    
        if(initCpus){
            parallel::stopCluster(cl)
        }

    }else{
        
        if(trace>0){
            requireNamespace("pbapply")
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

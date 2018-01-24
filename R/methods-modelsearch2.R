### methods-modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (16:43) 
## Version: 
## last-updated: jan 18 2018 (15:38) 
##           By: Brice Ozenne
##     Update #: 127
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * nStep
## ** documentation - nStep
#' @title Find the Number of Steps Performed During the Sequential Testing
#' @description Find the number of steps performed during the sequential testing
#' @name nStep
#' 
#' @param object a modelsearch2 object
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' nStep(res)
#'
#' @export
`nStep` <-
  function(object) UseMethod("nStep")

## ** function - nStep
#' @rdname nStep
#' @export
nStep.modelsearch2 <- function(object){

    return(length(object$sequenceTest))
    
}

## * getStep
## ** documentation - getStep
#' @title Extract one Step From the Sequential Procedure
#' @description Extract one step from the sequential procedure.
#' @name getStep
#' 
#' @param object a modelsearch2 object
#' @param step which test should be extracted?
#' @param slot the element from the modelsearch2 object that should be extracted.
#' @param ... not used.
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' getStep(res)
#' getStep(res, slot = "sequenceTest")
#' getStep(res, slot = "sequenceQuantile")
#' getStep(res, step = 1)
#'
#' @export
`getStep` <-
  function(object, ...) UseMethod("getStep")

## ** function - getStep
#' @rdname getStep
#' @export
getStep.modelsearch2 <- function(object, step = nStep(object), slot = NULL, ...){
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }
    if(!is.null(slot) && slot %in% names(object) == FALSE){
        stop("argument \'slot\' must be one of \"",paste(names(object),collapse="\" \""),"\"\n")
    }

    ## ** subset
    new.object <- object
    new.object$sequenceTest <- object$sequenceTest[step]
    new.object$sequenceModel <- object$sequenceModel[step]
    if(object$method.p.adjust == "max"){
        new.object$sequenceQuantile <- object$sequenceQuantile[[step]]
        sequenceIID <- object$sequenceIID[step]
        sequenceSigma <- object$sequenceSigma[step]
    }

    ## ** export
    if(is.null(slot)){
        return(new.object)
    }else{
        if(is.list(new.object[[slot]])){
            return(new.object[[slot]][[1]])
        }else{
            return(new.object[[slot]])
        }
    }
}

## * getNewLink
## ** documentation - getNewLink
#' @title Find the Links that Should be Added Accroding to the Sequential Testing
#' @description Find the links that should be added accroding to the sequential testing
#' @name getNewLink
#' 
#' @param object a modelsearch2 object
#' @param step which test should be extracted?
#' @param ... not used
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' getNewLink(res)
#'
#' @export
`getNewLink` <-
  function(object, ...) UseMethod("getNewLink")

## ** function - getNewLink
#' @rdname getNewLink
#' @export
getNewLink.modelsearch2 <- function(object, step = 1:nStep(object), ...){

    selected <- link <- NULL
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }

    ## ** extract
    ls.link <- lapply(step, function(x){
        iStep <- getStep(object,step=x,slot="sequenceTest")
        return(subset(iStep, subset = selected == TRUE, select = "link", drop = TRUE))
    })

    return(unlist(ls.link))    
}

## * merge
## ** documentation - merge
#' @title Merge two modelsearch Objects
#' @description Merge two modelsearch objects. Does not check for meaningful result.
#' @name merge
#' 
#' @param x a modelsearch2 object.
#' @param y a modelsearch2 object that will be added to x.
#' @param ... not used.
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res.x <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm", nStep = 2)
#' res.y <- modelsearch2(getStep(res.x, slot = "sequenceModel"), 
#'                       statistic = "score", method.p.adjust = "holm")
#' res.xy <- merge(res.x,res.y)
#'
#' modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")

## ** function - merge
#' @rdname merge
#' @export
merge.modelsearch2 <- function(x, y, ...){

    ## ** merge
    x$sequenceTest <- c(x$sequenceTest,y$sequenceTest)
    x$sequenceModel <- c(x$sequenceModel,y$sequenceModel)

    if(sum(c(y$method.p.adjust,x$method.p.adjust) == "max") == 1){
        if(x$method.p.adjust != "max"){
            lastStep.x <- nStep(x)
            Mtest.x <- getStep(x, step = 1, slot = "sequenceTest")
            nLink.x <- NROW(Mtest.x)
            name.x <- Mtest.x[["link"]]

            x$sequenceQuantile <- rep(NA, times = lastStep.x)
            if(!is.null(y$sequenceIID)){
                x$sequenceIID <- vector(mode = "list", length = lastStep.x)
            }
            x$sequenceSigma <- vector(mode = "list", length = lastStep.x)            
        }
        if(y$method.p.adjust != "max"){
            lastStep.y <- nStep(y)
            Mtest.y <- getStep(y, step = 1, slot = "sequenceTest")
            nLink.y <- NROW(Mtest.y)
            name.y <- Mtest.y[["link"]]

            y$sequenceQuantile <- rep(NA, times = lastStep.y)
            if(!is.null(x$sequenceIID)){
                y$sequenceIID <- vector(mode = "list", length = lastStep.y)
            }
            y$sequenceSigma <- vector(mode = "list", length = lastStep.y)            
        }
        x$sequenceQuantile <- c(x$sequenceQuantile,y$sequenceQuantile)
        x$sequenceIID <- c(x$sequenceIID,y$sequenceIID)
        x$sequenceSigma <- c(x$sequenceSigma,y$sequenceSigma)

    }
    
    for(iSlot in c("statistic","method.p.adjust","alpha","method.iid","cv")){
        x[[iSlot]] <- unique(c(x[[iSlot]],y[[iSlot]]))
        if(length(x[[iSlot]])>1){x[[iSlot]] <- NA}
    }
    ## ** export
    return(x)    
}


#----------------------------------------------------------------------
### methods-modelsearch2.R ends here

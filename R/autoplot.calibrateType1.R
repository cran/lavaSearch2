### autoplot.calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (13:20) 
## Version: 
## Last-Updated: apr  5 2018 (13:52) 
##           By: Brice Ozenne
##     Update #: 16
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
##' @title Graphical Display of the Bias or Type 1 Error
##' @description Graphical display of the bias or type 1 error
##' for the output of \code{\link{calibrateType1}}.
##' @name autoplot_calibrateType1
##' 
##' @param object output of \code{\link{calibrateType1}}.
##' @param type [character] if type equals \code{"bias"} the bias will be displayed.
##' Otherwise if it equals \code{"type1error"} the type 1 error will be displayed.
##' @param plot [logical] should the plot be displayed?
##' @param color.threshold [character] the color for the line representing the expected value(s).
##' @param type.bias [character] if type.bias equals \code{"absolute"} the absolute bias will be used.
##' Otherwise if it equals \code{"relative"} the relative bias will be used.
##' Only relevant when type equals \code{"bias"}.
##' @param alpha [numeric, 0-1] the significance threshold to consider.
##' Only relevant when type equals \code{"type1error"}.
##' @param nrow.legend [integer, >0] the number of rows for the legend.
##' Only relevant when type equals \code{"type1error"}.
##' @param name2label [named character vector] the label for the legend.
##' The vector should contain the method names (see details).
##' Only relevant when type equals \code{"type1error"}.
##' @param color [character vector] a vector of colours to be used to color the lines.
##' Only relevant when type equals \code{"type1error"}.
##' @param keep.method [character vector] the methods names for which the type 1 error should be displayed.
##' Only relevant when type equals \code{"type1error"}.
##'
##' @details Method names:
##' \itemize{
##' \item \code{p.Ztest}
##' \item \code{p.Satt}
##' \item \code{p.KR}
##' \item \code{p.robustZtest}
##' \item \code{p.robustSatt}
##' \item \code{p.robustKR}
##' }
##' @return An list containing:
##' \itemize{
##' \item plot: a ggplot object.
##' \item data: the dataset used to generate the ggplot object.
##' }
##' 


## * plotBias.calibrateType
##' @rdname autoplot_calibrateType1
##' @method autoplot calibrateType1
##' @export
autoplot.calibrateType1 <- function(object, type = "bias", plot = TRUE, color.threshold = "red",
                                    type.bias = "absolute",
                                    alpha = 0.05, nrow.legend = NULL, name2label = NULL, color = NULL, keep.method = NULL){

    type <- match.arg(type, choices = c("bias","type1error"))

    ## ** display bias
    if(type == "bias"){
        ## *** data
        type.bias <- match.arg(type.bias, choices = c("absolute","relative"))

        if(type.bias == "absolute"){
            object$bias$bias.ML <- object$bias$estimate.truth-object$bias$estimate.ML
            object$bias$bias.MLcorrect <- object$bias$estimate.truth-object$bias$estimate.MLcorrect
        }else if(type.bias == "relative"){
            object$bias$bias.ML <- object$bias$estiamte.truth-object$bias$estiamte.ML
            object$bias$bias.MLcorrect <- object$bias$estimate.truth-object$bias$estimate.MLcorrect
        }

        df1 <- cbind(object$bias[,c("n","rep","name","type","bias.ML")], corrected = FALSE)
        names(df1)[5] <- "bias"
        df2 <- cbind(object$bias[,c("n","rep","name","type","bias.MLcorrect")], corrected = TRUE)
        names(df2)[5] <- "bias"

        df.gg <- rbind(df1,df2)

        ## *** display
        gg <- ggplot() + geom_boxplot(data = df.gg, aes_string(x = "type", y = "bias")) 
        gg <- gg + facet_wrap(corrected ~ n, labeller = label_both)
        df.line <- data.frame(x = c(-Inf,Inf), y = 0)
        gg <- gg + geom_line(data = df.line, aes_string(x = "x", y = "y"), color = color.threshold)
        
    }

    ## ** display type 1 error
    if(type == "type1error"){
        ## *** data
        dfLong <- melt(object$p.value,
                       measure.vars = grep("^p.",names(object$p.value),value = TRUE),
                       value.name = "p.value",
                       variable.name = "method")
        df.gg <- stats::aggregate(dfLong$p.value,
                                  by = list(n = dfLong$n, method = dfLong$method, link = dfLong$link),
                                  FUN = function(x){c(n = length(x), type1error = mean(x<=alpha, na.rm = TRUE))},
                                  simplify = FALSE)
        df.gg <- cbind(df.gg[,c("n","method","link")],
                       do.call(rbind,df.gg[,"x"]))

        ## *** display
        if(is.null(keep.method)){
            keep.method <- as.character(unique(df.gg$method))
        }
        if(is.null(name2label)){
            name2label <- c(p.Ztest = "Gaussian approx.",
                            p.Satt = "Satterthwaite approx.",
                            p.KR = "Satterthwaite approx. with small sample correction",
                            p.robustZtest = "robust Gaussian approx.",
                            p.robustSatt = "robust Satterthwaite approx.",
                            p.robustKR = "robust Satterthwaite approx. with small sample correction"
                            )
        }
        if(is.null(color)){
            ## from ggthemes::colorblind_pal()(8)
            color <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
        }

        label.method <- name2label[keep.method]
        n.method <- length(keep.method)

        
        gg <- ggplot(df.gg, aes_string(x = "n", y = "type1error", group = "method", color = "method", shape = "method"))
        gg <- gg + geom_point(size = 3) + geom_line(size = 2)
        gg <- gg + facet_grid(~link, labeller = label_parsed)
        gg <- gg + geom_abline(intercept = alpha, slope = 0, color = color.threshold)
        gg <- gg + xlab("sample size")
        gg <- gg + ylab("type 1 error rate")
        gg <- gg + theme(legend.position = "bottom")

        if(!is.null(nrow.legend)){
            gg <- gg + guides(color=guide_legend(nrow=nrow.legend,byrow=TRUE))
        }
        
        gg <- gg + scale_color_manual("",
                                      breaks = keep.method,
                                      label = label.method,
                                      values = color[1:n.method])
        gg <- gg + scale_shape_manual("",
                                      breaks = keep.method,
                                      label = label.method,
                                      values = seq(15,by=1,length=n.method))
        
         
    }
    
    ## ** export
    if(plot){
        print(gg)
    }
    return(invisible(list(plot = gg,
                          data = df.gg)))
}


######################################################################
### autoplot.calibrateType1.R ends here

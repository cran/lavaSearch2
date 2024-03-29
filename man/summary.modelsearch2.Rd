% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelsearch2-summary.R
\name{summary.modelsearch2}
\alias{summary.modelsearch2}
\title{summary Method for modelsearch2 Objects}
\usage{
\method{summary}{modelsearch2}(object, print = TRUE, ...)
}
\arguments{
\item{object}{output of the \code{modelsearch2} function.}

\item{print}{should the summary be printed in the terminal.}

\item{...}{[internal] only used by the generic method.}
}
\description{
summary method for modelsearch2 objects.
}
\details{
The column \code{dp.Info} contains the percentage of extended models (i.e. model with one additional link)
for which the information matrix evaluated at the value of the parameters of the initial model is non positive definie.
}

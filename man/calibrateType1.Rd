% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrateType1.R
\name{calibrateType1}
\alias{calibrateType1}
\alias{calibrateType1.lvm}
\alias{calibrateType1.lvmfit}
\title{Simulation Study Assessing Bias and Type 1 Error}
\usage{
calibrateType1(object, param, n.rep, ...)

\method{calibrateType1}{lvm}(
  object,
  param,
  n.rep,
  n,
  correction = TRUE,
  warmup = NULL,
  null = NULL,
  F.test = FALSE,
  cluster = NULL,
  generative.object = NULL,
  generative.coef = NULL,
  true.coef = NULL,
  n.true = 1e+06,
  round.true = 2,
  bootstrap = FALSE,
  n.bootstrap = 1000,
  checkType1 = FALSE,
  checkType2 = FALSE,
  dir.save = NULL,
  label.file = NULL,
  seed = NULL,
  cpus = 1,
  trace = 2,
  ...
)

\method{calibrateType1}{lvmfit}(
  object,
  param,
  n.rep,
  correction = TRUE,
  F.test = FALSE,
  bootstrap = FALSE,
  n.bootstrap = 1000,
  seed = NULL,
  trace = 2,
  cpus = 1,
  ...
)
}
\arguments{
\item{object}{a \code{lvm} object defining the model to be fitted.}

\item{param}{[character vector] names of the coefficient whose value will be tested.}

\item{n.rep}{[integer, >0] number of simulations per sample size.}

\item{...}{[internal] Only used by the generic method.}

\item{n}{[integer vector, >0] sample size(s) considered in the simulation study.}

\item{correction}{[logical] should the type 1 error after correction be computed?}

\item{warmup}{[list of lvm] a list of \code{lvm} objects that will be sequentially fitted with for
starting values the parameter of the previous model in the list (if any). The parameters of the final
model of the list are used to initialize the fit of the model of interest (i.e. object).}

\item{null}{[numeric vector] vector of null hypotheses, one for each model coefficient.
By default a vector of 0.}

\item{F.test}{[logical] should a multivariate Wald test be perform testing simultaneously all the null hypotheses?}

\item{cluster}{[integer vector] the grouping variable relative to which the observations are iid.
Will be passed to \code{lava::estimate}.}

\item{generative.object}{[lvm] object defining the statistical model generating the data.}

\item{generative.coef}{[name numeric vector] values for the parameters of the generative model.
Can also be \code{NULL}: in such a case the coefficients are set to default values decided by lava (usually 0 or 1).}

\item{true.coef}{[name numeric vector] expected values for the parameters of the fitted model.}

\item{n.true}{[integer, >0] sample size at which the estimated coefficients will be a reliable approximation of the true coefficients.}

\item{round.true}{[integer, >0] the number of decimal places to be used for the true value of the coefficients. No rounding is done if \code{NULL}.}

\item{bootstrap}{[logical] should bootstrap resampling be performed?}

\item{n.bootstrap}{[integer, >0] the number of bootstrap sample to be used for each bootstrap.}

\item{checkType1}{[logical] returns an error if the coefficients associated to the null hypotheses do not equal 0.}

\item{checkType2}{[logical] returns an error if the coefficients associated to the null hypotheses equal 0.}

\item{dir.save}{[character] path to the directory were the results should be exported.
Can also be \code{NULL}: in such a case the results are not exported.}

\item{label.file}{[character] element to include in the file name.}

\item{seed}{[integer, >0] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
Can also be \code{NULL}: in such a case no seed is set.}

\item{cpus}{[integer >0] the number of processors to use.
If greater than 1, the simulations are performed in parallel.}

\item{trace}{[integer] should the execution of the function be trace. Can be 0, 1 or 2.}
}
\value{
An object of class \code{calibrateType1}.
}
\description{
Perform a simulation study over one or several sample size
to assess the bias of the estimate
and the type 1 error of the Wald test and robust Wald test
}
\examples{
\dontrun{
#### simulate data ####
m.Sim <- lvm(c(Y1[mu1:sigma]~1*eta,
               Y2[mu2:sigma]~1*eta,
               Y3[mu3:sigma]~1*eta,
               eta~beta1*Group+beta2*Gender))
latent(m.Sim) <- ~eta
categorical(m.Sim, labels = c("M","F")) <- ~Gender

d <- lava::sim(m.Sim, 1e2)

#### calibrate type 1 error on the estimated model ####
m <- lvm(Y1~eta,
         Y2~eta,
         Y3~eta,
         eta~Group+Gender)
e <- lava::estimate(m, data = d)
res <- calibrateType1(e, param = "eta~Group", n.rep = 100)
res <- calibrateType1(e, param = c("eta~Group","Y1~eta"), F.test = TRUE, n.rep = 100)
res <- calibrateType1(e, param = "eta~Group", n.rep = 100, cpus = 4)
summary(res)
}

}
\seealso{
\code{link{autoplot.calibrateType1}} for a graphical display of the bias or of the type 1 error.
}
\author{
Brice Ozenne
}

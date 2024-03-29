% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sCorrect-leverage2.R
\name{leverage2}
\alias{leverage2}
\alias{leverage2.lvmfit}
\alias{leverage2.lvmfit2}
\title{Leverage With Small Sample Correction.}
\usage{
leverage2(object, format, ssc, ...)

\method{leverage2}{lvmfit}(object, format = "wide", ssc = lava.options()$ssc, ...)

\method{leverage2}{lvmfit2}(object, format = "wide", ...)
}
\arguments{
\item{object}{a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).}

\item{format}{[character] Use \code{"wide"} to return the residuals in the wide format (one row relative to each sample).
Otherwise use \code{"long"} to return the residuals in the long format.}

\item{ssc}{[character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
Only relevant when using a \code{lvmfit} object.}

\item{...}{additional argument passed to \code{estimate2} when using a \code{lvmfit} object.}
}
\value{
a matrix containing the leverage relative to each sample (in rows)
and each endogenous variable (in column).
}
\description{
Extract leverage values from a latent variable model, with small sample correction.
}
\details{
The leverage are defined as the partial derivative of the fitted values with respect to the observations.
\deqn{
leverage_i = \frac{\partial \hat{Y}_i}{\partial Y_i}
}
See Wei et al. (1998). \cr \cr

When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the leverage.
}
\examples{
#### simulate data ####
set.seed(10)
m <- lvm(Y1~eta,Y2~eta,Y3~eta)
latent(m) <- ~eta
d <- lava::sim(m,20, latent = FALSE)

#### latent variable models ####
e.lvm <- estimate(m, data = d)
leverage2(e.lvm)

}
\references{
Bo-Cheng Wei et al., Generalized Leverage and its applications (1998), Scandinavian Journal of Statistics 25:1:25-37.
}
\seealso{
\code{\link{estimate2}} to obtain \code{lvmfit2} objects.
}
\concept{estimator}
\keyword{smallSampleCorrection}

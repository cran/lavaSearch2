% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntDensTri.R
\name{intDensTri}
\alias{intDensTri}
\title{Integrate a Gaussian/Student Density over a Triangle}
\usage{
intDensTri(
  mu,
  Sigma,
  df,
  n,
  x.min,
  z.max = NULL,
  type = "double",
  proba.min = 1e-06,
  prune = NULL,
  distribution = "pmvnorm"
)
}
\arguments{
\item{mu}{[numeric vector] the expectation of the joint distribution.}

\item{Sigma}{[matrix] the variance-covariance of the joint distribution.}

\item{df}{[integer > 0] the degree of freedom of the joint Student's t distribution.
Only used when \code{distribution="pvmt"}.}

\item{n}{[integer > 0] number of points for the numerical integration.}

\item{x.min}{[numeric] the minimum value along the x axis.}

\item{z.max}{[numeric vector, optional] the maximum value along the z axis.
Define the dimension of Z.}

\item{type}{[character] the type of mesh to be used.
Can be \code{\"raw\"}, \code{\"double\"}, or \code{\"fine\"}.}

\item{proba.min}{[numeric 0-1] the probability used to find the maximum value along the x axis.
Only used if \code{prune} is not specified.}

\item{prune}{[integer >0] number of standard deviations after which the domain ends along the x axis.}

\item{distribution}{[character] type of joint distribution.
Can be \code{"pmvnorm"} (normal distribution) or \code{"pvmt"} (Student's t distribution)}
}
\value{
A numeric.
}
\description{
Consider a univariate random variable X,
two multivariate random variables Y and Z,
and t1 and t2 two real numbers.
This function can compute either
P[|X|>t1,|X]>|Y1|,...,|X]>|Yp|] if zmin is not specified,
P[|Z1|<t2,...,|Zq|<t2,|X|>t1,|X]>|Y1|,...,|X]>|Yp|] if zmin is specified.
}
\details{
Argument \code{type}: \itemize{
\item \code{\"raw\"}: mesh with points inside the domain
\item \code{\"double\"}: mesh with points outside the domain
\item \code{\"fine\"}: mesh with points inside the domain plus additional rectangles trying to fill the missing domain.
}

Argument \code{Sigma} and \code{mu}:
define the mean and variance-covariance of the random variables X, Y, Z
(in this order). The length of the argument \code{z.max} is used to define the dimension of Z.
The dimension of X is always 1.
}
\examples{
library(mvtnorm)

p <- 2
Sigma <- diag(p)
mu <- rep(0, p)

## bivariate normal distribution
z2 <- qmvt(0.975, mean = mu, sigma = Sigma, df = 1e3)$quantile

# compute integral
intDensTri(mu = mu, Sigma = Sigma, n=5, x.min=0, type = "fine")$value-1/2
intDensTri(mu = mu, Sigma = Sigma, n=30, x.min=0, type = "raw")$value-1/2
intDensTri(mu = mu, Sigma = Sigma, n=50, x.min=0, type = "raw")$value-1/2

intDensTri(mu = mu, Sigma = Sigma, df = 5, n=5, x.min=0, distribution = "pmvt")$value-1/2
res <- intDensTri(mu = mu, Sigma = Sigma, df = 5, n=10, x.min=0, distribution = "pmvt")
res$value-1/2
ggplot2::autoplot(res)

## trivariate normal distribution
\dontrun{
p <- 3
Sigma <- diag(p)
mu <- rep(0, p)

res2 <- intDensTri(mu = mu, Sigma = Sigma, n=5, x.min = 0, z.max = 10)
ggplot2::autoplot(res2)
ggplot2::autoplot(res2, coord.plot = c("x","z1"))
res2
}

#### when the distribution is far from 0
\dontrun{
eq1 <- intDensTri(mu = c(10,0), Sigma = diag(1,2), 
                  x.min = 2, n=10)
eq1$value-1
ggplot2::autoplot(eq1)

eq2 <- intDensTri(mu = c(10,0,0), Sigma = diag(1,3),
                  x.min=2, z.max = 10, type = "raw",
                  n=10)
ggplot2::autoplot(eq2, coord.plot = c("y1","z1"))
eq2$value-1

## more variables
p <- 5
Sigma <- diag(p)
mu <- rep(0, p)

res2 <- intDensTri(mu = mu, Sigma = Sigma, n=5, x.min = 1, z.max = c(2,2))
res2$grid
}

}
\concept{post-selection inference}

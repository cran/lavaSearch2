% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information2.R
\name{information2-internal}
\alias{information2-internal}
\alias{.information2}
\alias{.hessian2}
\title{Compute the Expected Information Matrix From the Conditional Moments}
\usage{
.information2(
  dmu,
  dOmega,
  Omega,
  n.corrected,
  index.Omega,
  leverage,
  n.cluster,
  grid.meanparam,
  n.grid.meanparam,
  grid.varparam,
  n.grid.varparam,
  name.param,
  n.param
)

.hessian2(
  dmu,
  d2mu,
  dOmega,
  d2Omega,
  Omega,
  n.corrected,
  index.Omega,
  leverage,
  n.cluster,
  grid.meanparam,
  n.grid.meanparam,
  grid.varparam,
  n.grid.varparam,
  name.param,
  n.param,
  residuals
)
}
\description{
Compute the expected information matrix from the conditional moments.

Compute the Hessian matrix from the conditional moments.
}
\details{
\code{.information2} will perform the computation individually when the
argument \code{index.Omega} is not null.

\code{.hessian} will perform the computation individually when the
argument \code{index.Omega} is not null.
}
\keyword{internal}

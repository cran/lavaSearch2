% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information2.R
\name{dInformation2-internal}
\alias{dInformation2-internal}
\alias{.dInformation2}
\title{Compute the First Derivative of the Expected Information Matrix}
\usage{
.dInformation2(
  dmu,
  d2mu,
  dOmega,
  d2Omega,
  Omega,
  OmegaM1,
  n.corrected,
  index.Omega,
  leverage,
  n.cluster,
  name.param,
  name.3deriv
)
}
\description{
Compute the first derivative of the expected information matrix.
}
\details{
\code{.dInformation2} will perform the computation individually when the
argument \code{index.Omega} is not null.
}
\keyword{internal}

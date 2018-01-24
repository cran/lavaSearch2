#' @docType package
#' @name lavaSearch
#'
#' @title Tools for Model Specification in the Latent Variable Framework
#' @description Tools for model specification in the latent variable framework
#' (add-on to the lava package). The package contains three main functionalites:
#' Wald tests/F-tests with improved control of the type 1 error in small samples,
#' adjustment for multiple comparisons when searching for local dependencies,
#' and adjustment for multiple comparisons when doing inference for multiple latent variable models. 
#'
#' @details The latent variable models (lvm) considered in this package can be written \cr
#' as a measurement model:
#' \deqn{Y_i = \nu + \eta_i \Lambda + X_i K + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + \eta_i B + X_i \Gamma + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr \cr
#' 
#' The corresponding conditional mean is:
#' \deqn{
#' \mu_i(\theta) = E[Y_i|X_i] = \nu + (\alpha + X_i \Gamma) (1-B)^{-1} \Lambda + X_i K
#' }
#' \deqn{
#' \Omega(\theta) = Var[Y_i|X_i] = \Lambda^{t} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda + \Sigma
#' }
#'
#' Therefore:
#' \itemize{
#' \item \eqn{\nu}, \eqn{K}, \eqn{\alpha}, \eqn{\Gamma} are pure mean parameters.
#' \item \eqn{\Psi}, \eqn{\Sigma} pure variance parameters.
#' \item \eqn{\Lambda}, \eqn{B} are both mean and variance parameters.
#' }
#' 

#' @import lava
#' @import ggplot2
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom graphics par plot text
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag
#' @importFrom methods is
#' @importFrom multcomp glht 
#' @importFrom mvtnorm pmvnorm qmvnorm rmvnorm qmvt pmvt
#' @importFrom numDeriv jacobian hessian
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom reshape2 melt
#' @importFrom stats anova as.formula coef cov df.residual dist formula hclust logLik median model.frame model.matrix na.omit optim p.adjust pf pnorm predict qqnorm quantile pt residuals rnorm sd setNames sigma update vcov
#' @importFrom tcltk setTkProgressBar tkProgressBar
#' @importFrom utils methods packageVersion setTxtProgressBar tail txtProgressBar
#' 
NULL



  

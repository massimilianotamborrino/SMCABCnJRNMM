#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL


#'@rdname model
#'@title model
#'@description (faster version of the) Strang splitting method for path simulation of the
#' observed output process from the stochastic multi-population JRNMM. This is analogous
#' to the observedJRNMM routine in the StrangSplitting R package
#'@param N number of populations of neurons
#'@param grid time points at which to simulate the process
#'@param h step size used for path simulation
#'@param startv starting point x0
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param dSigma vector of diagonal entries of the Sigma matrix
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'#'@return Observed output process from the stochastic N-population JRNMM
#'@export
model <- function(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K){
  Y<-fast_JRNMM_Splitting(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K)
  return(model_(N,Y))
}

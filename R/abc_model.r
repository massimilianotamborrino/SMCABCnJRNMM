#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL


#'@rdname model
#'@title model
#'@description (faster version of the) Strang splitting method for path simulation of the
#' stochastic multi-population JRNMM
#'@param N number of populations of neurons
#'@param grid time points at which to simulate the process
#'@param h step size used for path simulation
#'@param startv starting point x0
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param dSigma vector of diagonal entries of the Sigma matrix
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'#'@return path of the stochastic N-population JRNMM
#'@export
model <- function(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K){
  Y<-fast_JRNMM_Splitting(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K)
  return(model_(N,Y))
}

#'@rdname KmatrixgivenLc
#'@title Code to calculate the Kmatrix given L and c
#'@description Sampling from the Gaussian Perturbation kernel for the continuous parameters
#'@param N number of populations of neurons
#'@param L strength parameter
#'@param c [0,1] parameter
#'@return Kmatrix given L and c
#'@export
KmatrixgivenLc <- function(N,L,c){
  return(KmatrixgivenLc_(N,L,c))
}

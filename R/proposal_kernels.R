#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL


#'@rdname perturb_continuous
#'@title perturb_continuous
#'@description Sampling from the Gaussian Perturbation kernel for the continuous parameters
#'@param Theta  Mean of the Gaussian perturbation kernel
#'@param Sigma covariance matrix of the Gaussian perturbation kernel
#'@return 1 sampled value from the Gaussian Perturbation kernel for the continuous parameters
#'@export
perturb_continuous <- function(Theta,Sigma){
  return(perturb_continuous_(Theta,Sigma))
}

#'@rdname perturb_discrete
#'@title perturb_discrete
#'@description Sampling from a Bernoulli distribution, with a further perturbation on the sampled value according to a q_stay probability
#'@param N number of neural populations
#'@param pvec N*(N-1)-dimensional probability vector of the binary parameter
#'@param stay_prob Probability of not perturbing the sampled binary value
#'@return 1 sampled value from the discrete Perturbation kernel for the discrete parameters
#'@export
perturb_discrete <- function(N,pvec,stay_prob){
  Rho<-matrix(Inf, nrow=N,ncol=N)
  counter<-0
  outcomes<-c(1,0)
  pert<-c(0,1)
  sam<-rep(0,N*(N-1))
  for(i in 1:(N*(N-1))) sam[i]<- sample(outcomes,1,prob=c(pvec[i],1-pvec[i]))
  samp.pert<-sample(pert,N*(N-1),replace=T,prob=c(stay_prob,1-stay_prob))
  res<-abs(sam-samp.pert)
  for(j in 1:N){
    for(k in 1:N){
      if(k!=j){
        counter<-counter+1
        Rho[j,k]<-res[counter]
      }
    }
  }
  return(Rho)}


 

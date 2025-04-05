#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp evalCpp
NULL

#'@rdname problemprior
#'@title problemprior
#'@description Either pdf of priors (for draw=0) or sample from the prior (for draw =1)
#'@param theta Parameter Vector;
#'@param draw Draw=0 if we want to compute the product of the pdfs of the priors
#'Draw =1 if we want to sample from the priors
#'@param Pr_cont Matrix of lower/upper bounds for the uniform priors for the continuous parameters
#'@param N  number of populations of neurons
#'@param whichprior 1 uniform
#'@return Simulation of observations from all models (for draw=1) or pdf of the priors (for draw=0), returned as vector;
#'@export
problemprior<-function(theta,draw,Pr_cont, N,whichprior='unif') {
  {if(whichprior=='unif') {return(nJRNMM_prior_(theta,draw,Pr_cont,N));}
   # else if(whichprior=='lognormal') {return(FHN_prior2_(theta,draw));}
  #  else return(FHN_prior3_(theta,draw))}
  }
}


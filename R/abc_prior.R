#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp evalCpp
NULL

#'@rdname problemprior
#'@title problemprior
#'@description Either pdf of priors (for draw=0) or sample from the prior (for draw =1)
#'@param theta Parameter Vector;
#'@param draw Draw=0 if we want to compute the product of the pdfs of the priors
#'Draw =1 if we want to sample from the priors
#'@param whichprior 1 uniform, 2 lognormal, 3 exponential
#'@param Pr_cont Matrix of lower/upper bounds for the uniform priors for the continuous parameters
#'@return Simulation of observations from all models (for draw=1) or pdf of the priors (for draw=0), returned as vector;
#'@export
problemprior<-function(theta,draw,Pr_cont,whichprior='unif') {
  {if(whichprior=='unif') {return(nJRNMM_prior_(theta,draw,Pr_cont));}
   # else if(whichprior=='lognormal') {return(FHN_prior2_(theta,draw));}
  #  else return(FHN_prior3_(theta,draw))}
  }
}


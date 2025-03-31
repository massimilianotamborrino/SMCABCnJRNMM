#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname abc_distance
#'@title abc_distance
#'@description Compute the distance given the summaries
#'@param summaries_parameters vector of summaries
#'@param summaries_weights vector with weights of summaries computed from the data
#'@param N  number of populations of JRNMMs
#'@return ABC distance for the given summaries
#'@export
abc_distance<-function(summaries_parameters,summaries_weights,N) {

  #access summary weights
  #we_dens<-summaries_weights[1]
 # we_coh<-summaries_weights[2] REMOVED
#we_ccf<-summaries_weights[3]


  return((summaries_parameters[1]*summaries_weights[1]+summaries_parameters[2]*summaries_weights[2]+#summaries_parameters[3]/((N-1)/2)*summaries_weights[2]+
            summaries_parameters[3]/(N-1)*summaries_weights[3])/N)
}



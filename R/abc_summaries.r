#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL

#-------------------------------------------------------------------------
#'@rdname abc_summaries
#'@title abc_summaries
#'@description Computation of the three ABC summary statistics
#'@param N number of populations of JRNMMs
#'@param summaries_X            summaries of observed dataset X
#'@param summaries_Y            summaries of synthetic dataset X
#'@param summaries_extra  vector of quantities needed to compute the summary statistics
#'@return distance value among two sets of summaries
#'@export
abc_summaries<-function(N,summaries_X,summaries_Y,summaries_extra)
{

  #access the different summaries
  #densities
  dens_matX<-summaries_X[[1]]
  dens_matY<-summaries_Y[[1]]
  #spectral densities
  spec_matX<-summaries_X[[2]]
  spec_matY<-summaries_Y[[2]]
  ##coh's
  #coh_matX<-summaries_X[[3]]
  #coh_matY<-summaries_Y[[3]]
  #ccf's
  ccf_matX<-summaries_X[[3]]
  ccf_matY<-summaries_Y[[3]]

  #access specific summary parameters
  stepD<-summaries_extra[4]
  stepP<-summaries_extra[7]
  stepCcf<-summaries_extra[10]

  #Compute distance
  dist_dens<-stepD*sum(abs(dens_matX-dens_matY))
  dist_spec<-stepP*sum(abs(spec_matX-spec_matY))
  # dist_coh<-stepP*sum(abs(coh_matX-coh_matY))
  dist_ccf<-stepCcf*sum(abs(ccf_matX-ccf_matY))

  return(c(dist_dens,dist_spec,dist_ccf))
#  return(c(dist_dens,dist_spec,dist_coh,dist_ccf))
}


#' #-------------------------------------------------------------------------
#'@rdname abc_summaries_extra
#'@title abc_summaries_extra
#'@description Determination of quantities needed to determine ABC summaries of a given dataset
#'@param X       Nxn matrix of observed dataset/path of N populations of JR-NMMs with n points each
#'@param T       time horizon for path simulation
#'@param h       step size for path simulation
#'@return vector of required quantities for then calculating the summary statistics
#'@export
abc_summaries_extra<-function(X,T,h){

  #parameters for density
  startSupp<--40
  endSupp<-40
  Lsupport<-1001
  stepD<-(endSupp-startSupp)/(Lsupport-1)

  #parameters for spectral densities and MSCs
  span_val<-5*T
  lag_val<-100

  #X12<-ts.union(ts(X[1,],frequency = 1/h),ts(X[2,],frequency = 1/h))

  specX<-spectrum(X[1,],log="no",span=span_val,plot=FALSE) #specX<-spectrum(X12,spans = span_val)#crossSpectrum(X12,spans = span_val)
  spx<-specX$freq/h
  stepP<-diff(spx)[1]
  Lspec<-length(spx)

  #parameters for cross-correlation function
  Lccf<-201
  stepCcf<-h

  #compute the vector of summary parameters and return it
  summaries_extra<-c(startSupp,endSupp,Lsupport,stepD,span_val,Lspec,stepP,lag_val,Lccf,stepCcf)
  return(summaries_extra)
}

#-------------------------------------------------------------------------
#'@rdname abc_summaries_X
#'@title abc_summaries_X
#'@description Calculations of the ABC summaries of the reference data X
#'@param X       Nxn matrix of observed dataset/path of N populations of JR-NMMs with n points each
#'@param T       time horizon for path simulation
#'@param h       step size for path simulation
#'@param summaries_extra vector of quantities needed to compute the summary statistics
#'@return summaries (densities, spectral densities, MSCs, ccfs) of a dataset X
#'@export
abc_summaries_X<-function(X,T,h,summaries_extra)
{

  N<-dim(X)[1] #number of populations in the JRNMM

  #access specific summary parameters
  startSupp<-summaries_extra[1]
  endSupp<-summaries_extra[2]
  Lsupport<-summaries_extra[3]
  span_val<-summaries_extra[5]
  Lspec<-summaries_extra[6]
  lag_val<-summaries_extra[8]
  Lccf<-summaries_extra[9]
  #determine densties, spectral densities, coh's (MSCs) and ccf's of X
  dens_mat<-matrix(0,nrow=N,ncol=Lsupport)
  spec_mat<-matrix(0,nrow=N,ncol=Lspec)
  # coh_mat<-matrix(0,nrow=(N*(N-1)/2),ncol=Lspec)
  ccf_mat<-matrix(0,nrow=(N*(N-1)),ncol=Lccf)

  # count_coh<-0
  count_ccf<-0
  for(j in 1:(N-1)){
    Xj<-X[j,]
    spec_mat[j,]<-spectrum(Xj,log="no",span=span_val,plot=FALSE)$spec*h
    dens_mat[j,]<-density(Xj,n=Lsupport,from=startSupp,to=endSupp)$y
    for(k in (j+1):N){
      # Xjk<-ts.union(ts(X[j,], frequency = 1/h),ts(X[k,], frequency = 1/h))
      #  specXjk<- crossSpectrum(Xjk,spans = span_val)

      #coh
      #  count_coh<-count_coh+1
      # coh_mat[count_coh,]<-specXjk$coh

      #ccf
      Xk<-X[k,]
      count_ccf<-count_ccf+1
      ccfXjk<-ccf(Xj,Xk,lag.max=lag_val,plot=FALSE)
      ccf_mat[count_ccf,]<-ccfXjk$acf
      count_ccf<-count_ccf+1
      ccfXkj<-ccf(Xk,Xj,lag.max=lag_val,plot=FALSE)
      ccf_mat[count_ccf,]<-ccfXkj$acf
    }
    #spec
    # spec_mat[j,]<-specXjk$spec1
    #dens
  }
  spec_mat[N,]<-spectrum(Xk,log="no",span=span_val,plot=FALSE)$spec*h
  dens_mat[N,]<-density(Xk,n=Lsupport,from=startSupp,to=endSupp)$y
  #spec_mat[N,]<-specXjk$spec2
  #dens_mat[N,]<-density(X[N,],n=Lsupport,from=startSupp,to=endSupp)$y

  ret<-list()
  ret[[1]]<-dens_mat
  ret[[2]]<-spec_mat
  #  ret[[3]]<-coh_mat
  ret[[3]]<-ccf_mat

  #return summaries
  return(ret)
}


#-------------------------------------------------------------------------
#'@rdname abc_summaries_weights
#'@title abc_summaries_weights
#'@description Computation of weights (for distance measure among two sets of summaries) from observed data
#'@param summaries  summary statistics of observed vs simulated dataset
#'@param summaries_extra vector of quantities needed to compute the summary statistics
#'@param N number of populations in the JRNMM
#'@return  vector of summary weights
#'@export
abc_summaries_weights<-function(summaries,summaries_extra,N){

  #access specific summary parameters
  stepP<-summaries_extra[7]
  stepCcf<-summaries_extra[10]

  #access different summaries
  dens_mat<-summaries[[1]] #dens
  spec_mat<-summaries[[2]] #spec
  # coh_mat<-summaries[[3]] #coh
  ccf_mat<-summaries[[3]] #ccf

  specA<-stepP*sapply(1:N,FUN=function(I) sum(abs(spec_mat[I,])))
  #  cohA<-stepP*sapply(1:((N*(N-1))/2),FUN=function(I) sum(abs(coh_mat[I,])))
  ccfA<-stepCcf*sapply(1:(N*(N-1)),FUN=function(I) sum(abs(ccf_mat[I,])))

  mean_specA<-mean(specA)
  # mean_cohA<-mean(cohA)
  mean_ccfA<-mean(ccfA)

  #weight density
  we_dens<-mean_specA

  #  #weight coh
  # we_coh<-mean_specA/mean_cohA

  #weight ccf
  we_ccf<-mean_specA/mean_ccfA

  #determine vector of summaries and return it
  #summaries_weights<-c(we_dens,we_coh,we_ccf)
  summaries_weights<-c(we_dens,1,we_ccf)
  return(summaries_weights)
}



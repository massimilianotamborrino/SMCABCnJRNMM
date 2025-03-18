#' @useDynLib SMCABCnJRNMM
#' @importFrom Rcpp sourceCpp
NULL


#'@rdname SMCABCnJRNMM_acceptrate
#'@title SMCABCnJRNMM_acceptrate
#'@description SMC-ABC for time-varying thresholds computed as percentiles of ACCEPTED distances with a fixed budget (number of simulations)
#'@param data observed data
#'@param Pr_cont Matrix of lower/upper bounds for the uniform priors for the continuous parameters
#'@param stay_prob Probability of not perturbing the sampled discrete connection probability
#'@param extra_model List with Number of points to simulate N, final time T, grid, time step h, initial position, vector of diagonal entries of the gamma matrix, vector of diagonal entries of the sigma matrix, Theta0
#'@param ABCthreshold Initial value for the threshold for the guided approach
#'@param minaccept.rate Lowest acceptance rate after which we stop the algorithm
## #'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param numparticles number of particles to keep
#'@param alpha percentile to compute the ABC threshold out of all distances
#'@param sampling specify the chosen samplers. The options are standard or olcm
#'@param attempt number of iteration
#'@param folder specify the folder where the results are going to be saved in
#'@param whichprior choose between 'unif','lognormal' and and 'exp'
#'@param subsamplingby every how many simulated points the observation should be taken. The default=1, i.e., no subsampling
#'@export
#'
SMCABCnJRNMM_acceptrate<- function (data, Pr_cont,stay_prob,extra_model, ABCthreshold, minaccept.rate, numparticles, alpha, sampling,
                             attempt, folder, whichprior = 'unif',subsamplingby=1)
{
  N<-extra_model[[1]]
  T<-extra_model[[2]]
  grid<-extra_model[[3]]
  h<-extra_model[[4]]
  startv<-extra_model[[5]]
  dGamma<-extra_model[[6]]
  dSigma<-extra_model[[7]]
  Theta0<-extra_model[[8]]
  ndata<-length(grid)

  #% THE OBSERVED SUMMARIES - Code for the FHN
  summaries_extra<-abc_summaries_extra(data,T,h)
  summobs<-abc_summaries_X(data,T,h,summaries_extra)
  summ_weights<-abc_summaries_weights(summobs,summaries_extra,N)
  #####
  nfreepar_c <- N+2# % number of continuous parameters
  nfreepar_d <- N^2# % number of discrete parameters
  nfreepar <- nfreepar_c+nfreepar_d# % total number of parameters to be inferred (-N in fact)

  ABCdraws <- matrix(0,nrow=nfreepar,ncol=numparticles); #% will store accepted parameters (for 1 iteration only, not all)
  nsummaries <- length(summ_weights)#% the number of summary statistics
  distance_accepted <- matrix(0,nrow=1,ncol=numparticles); #% this is only for the olcm proposal
  lengthmodel<-length(data)
  #% initialization: t is the iteration counter
  t <- 1;
  numproposals0<-0
  number_sim<- numparticles/minaccept.rate
  tic()   #% This will  count the seconds required to obtain the desired number of accepted particles for each iteration
  RES<-foreach(success = 1:numparticles,.combine='rbind',.packages=c('SMCABCnJRNMM')) %dopar% {
    distance<-ABCthreshold+1
    numproposals <- 0;
    #simsumm_all <- c();

    while(distance >= ABCthreshold & numproposals<number_sim){
      numproposals <- numproposals +1;
      parameters <-  problemprior(-1,1,Pr_cont,whichprior);# % propose from the prior
      Theta0[,1]<- parameters[1:N]
      K<-KmatrixgivenLc(N,parameters[N+1],parameters[N+2])
      Rho<-matrix(parameters[-(1:nfreepar_c)],nrow=N,byrow=T)
      simdata <- model(N, grid, h, startv, dGamma,dSigma, Theta0,Rho,K); #// simulate from the model
      simdata<- simdata[,seq(1,ndata,by=subsamplingby)]
      summaries_Y<-abc_summaries_X(simdata,T,h,summaries_extra)
      simsumm <- abc_summaries(N,summobs,summaries_Y,summaries_extra)#
      distance <- abc_distance(simsumm,summ_weights,N)#
     # simsumm_all <-  cbind(simsumm_all,simsumm)
    }
    if(numproposals>=number_sim) {list(numproposals);stop('a particle got stucked');}
    return(list(numproposals,distance,parameters))#    return(list(numproposals,distance,parameters,simsumm_all))
  }
  eval_time <- toc()
  weights <- matrix(1,ncol=numparticles);
  norm_weights_c <- weights/sum(weights); # % normalised weights
  ess <- 1/sum(norm_weights_c^2); # % the Effective Sample Size

  numproposals<-sum(unlist(RES[,1]))
  distance_accepted<-matrix(as.numeric(unlist(RES[,2])),nrow=1,ncol=numparticles)
  ABCdraws<-matrix(as.numeric(unlist(RES[,3])),nrow=nfreepar,ncol=numparticles)
  ABCdraws_c<-ABCdraws[1:nfreepar_c,]
  ABCdraws_d<-ABCdraws[(nfreepar_c+1):nfreepar,]
#  simsumm_all<-matrix(as.numeric(unlist(RES[,4])),nrow=nsummaries)
  rm(RES)

  write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ABCthreshold,file=sprintf('%s/ABCthreshold_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ABCdraws_c,file=sprintf('%s/ABCdraws_c_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ABCdraws_d,file=sprintf('%s/ABCdraws_d_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(norm_weights_c,file=sprintf('%s/normweights_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)


  ABCthreshold_new<- ABCthreshold
  acceptrate<-numparticles/numproposals
  hat_p_vec<-sapply(1:N^2,FUN=function(I) sum(ABCdraws_d[I,])/numparticles)

  while (acceptrate>minaccept.rate) {#(totnumproposals < number_sim){ # I stop when the acceptance rate goes below a given minimum amount
    t <- t+1
    {if(acceptrate< minaccept.rate*10){ABCthreshold_temp <- quantile(distance_accepted,alpha/100);
    }
  else ABCthreshold_temp <- median(distance_accepted);}
    ABCthreshold_old <- ABCthreshold_new;
    ABCthreshold_new <- min(ABCthreshold_temp,ABCthreshold_old)#
    if (ABCthreshold_new == ABCthreshold_old)  {ABCthreshold_new <- 0.95*ABCthreshold_new;
    sprintf('Forced decrease of ABC threshold to: %f',ABCthreshold)}


    if(sampling=="standard"){
      C <- t(ABCdraws_c) - matrix(rep(norm_weights_c %*% t(ABCdraws_c), numparticles), nrow=numparticles,byrow=T);#   % subtract weighted mean                                                  % Remove mean (which is, also, weighted)
      C <- t(C) %*% (C * matrix(rep(t(norm_weights_c),nfreepar_c),nrow=numparticles));       # % Weighted Covariance Matrix
      C <- C / (1-sum(norm_weights_c^2)); #% the weighted covariance matrix
      C <- 0.5 * (C + t(C));   #% ensure symmetry
      Sigma <- 2*C;
      write.table(Sigma,file=sprintf('%s/Sigma_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)

    }
    else if (sampling=="olcm"){
      cov_olcm <- matrix(0,nrow=nfreepar_c,ncol=nfreepar_c);
      id_olcm <- which(distance_accepted < ABCthreshold_new);
      N0 <- length(id_olcm);  #% number of particles that satisfy the newest threshold
      if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
      { print('There are no particles satisfying the OLCM criterion. We must leave..')
        return(3)}
      weights_olcm <- weights[id_olcm];
      normweights_olcm <- weights_olcm/sum(weights_olcm);
      cov_olcm_all <- c();
      for (ii in 1:numparticles)
      {
        for (jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws_c[,id_olcm[jj]]-ABCdraws_c[,ii])%*%t(ABCdraws_c[,id_olcm[jj]]-ABCdraws_c[,ii]);
        cov_olcm <- (cov_olcm+t(cov_olcm))/2;
        if(isposdef(cov_olcm)==0) cov_olcm<-nearPD(cov_olcm,base.matrix=TRUE)$mat
        cov_olcm_all <- cbind(cov_olcm_all,cov_olcm);
      }
      write.table(cov_olcm_all,file=sprintf('%s/Sigma_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    }

    tic()

    RES2<-foreach(success = 1:numparticles,.combine='rbind',.packages=c('SMCABCnJRNMM')) %dopar% {
      distance<- ABCthreshold_new+1
      numproposals <- 0;
      numproposals0 <- 0;
      numproposalsneg <- 0;

#   simsumm_all<- c();
      while(distance>=ABCthreshold_new & numproposals<number_sim){

        index <- sample(1:length(norm_weights_c),size=1,prob=norm_weights_c);#needed to compute the mean of the proposal sampler

        if(sampling=='standard') theta_c <- perturb_continuous(ABCdraws_c[,index],Sigma)
        else if(sampling=='olcm') {
          cov_olcm <-matrix(0, nrow=nfreepar_c,ncol=nfreepar_c);
          for(jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws_c[,id_olcm[jj]]-ABCdraws_c[,index])%*%t(ABCdraws_c[,id_olcm[jj]]-ABCdraws_c[,index]);
          cov_olcm <- (cov_olcm+t(cov_olcm))/2;
          if(isposdef(cov_olcm)==0) cov_olcm<-nearPD(cov_olcm,base.matrix=TRUE)$mat  # This IF statement is useful if the proposal_cov is not definite positive
          #    % the above covariance is not "global" but is instead specific for the sampled particle
          theta_c <- perturb_continuous(ABCdraws_c[,index],cov_olcm)}
        if(min(theta_c)<0) numproposalsneg<-numproposalsneg+1
        prior <- problemprior(theta_c,0,Pr_cont,whichprior); #% evaluate prior
        if(prior==0) numproposals0<-numproposals0+1
        else {
          numproposals <- numproposals +1;
          #sample from the discrete kernel and perturb it according to q_stay
          Theta0[,1]<-theta_c[1:N]
          K<-KmatrixgivenLc(N,theta_c[N+1],theta_c[N+2])
          Rhov<-perturb_sample_discrete_(N,hat_p_vec,stay_prob)
          Rhom<-matrix(Rhov,nrow=N,byrow = T)
          simdata <- model(N, grid, h, startv, dGamma,dSigma, Theta0,Rhom,K); #// simulate from the model
          simdata<- simdata[,seq(1,ndata,by=subsamplingby)]
          summaries_Y<-abc_summaries_X(simdata,T,h,summaries_extra)
          simsumm <- abc_summaries(N,summobs,summaries_Y,summaries_extra)#
          distance <- abc_distance(simsumm,summ_weights,N)#
      #    simsumm_all <-  cbind(simsumm_all,simsumm);
        }#
      }
      if(numproposals>=number_sim) {list(numproposals);stop('a particle got stucked');}
      if(sampling=='standard'){
        dens <-0;
        for (ii in 1:numparticles) dens <- dens + norm_weights_c[ii]*dmvn(theta_c,ABCdraws_c[,ii],Sigma);
      }
      else if(sampling=='olcm'){
        dens <-0;
        for (ii in 1:numparticles) dens <- dens + norm_weights_c[ii]*dmvn(theta_c,ABCdraws_c[,ii],cov_olcm_all[,((ii-1)*nfreepar_c+1):(ii*nfreepar_c)]);
      }
     if(dens==0) return(print('error'))
      list(numproposals,distance,theta_c,Rhov,prior/dens,numproposals0,numproposalsneg)#  list(numproposals,distance,theta_c,Rhov,simsumm_all,prior/dens,numproposals0,numproposalsneg)
    }
    eval_time <- toc()

    numproposals<-sum(unlist(RES2[,1]))
    if(numproposals>number_sim) return(numproposals)
    distance_accepted<-matrix(as.numeric(unlist(RES2[,2])),nrow=1,ncol=numparticles)
    ABCdraws_c<-matrix(as.numeric(unlist(RES2[,3])),nrow=nfreepar_c,ncol=numparticles)
    ABCdraws_d<-matrix(as.numeric(unlist(RES2[,4])),nrow=nfreepar_d,ncol=numparticles)
  #  simsumm_all<-matrix(as.numeric(unlist(RES2[,5])),nrow=nsummaries)
    weights<-matrix(as.numeric(unlist(RES2[,5])),ncol=numparticles);
    numproposals0<-sum(unlist(RES2[,6]))
    numproposalsneg<-sum(unlist(RES2[,7]))
    rm(RES2)

#    totnumproposals<- totnumproposals+numproposals

    norm_weights_c <- weights/sum(weights); # % normalised weights
    ess <- 1/sum(norm_weights_c^2); # % the Effective Sample Size

    write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ABCthreshold_new,file=sprintf('%s/ABCthreshold_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ABCdraws_c,file=sprintf('%s/ABCdraws_c_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ABCdraws_d,file=sprintf('%s/ABCdraws_d_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposals0,file=sprintf('%s/numproposals0_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposalsneg,file=sprintf('%s/numproposalsneg_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(norm_weights_c,file=sprintf('%s/normweights_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)

    acceptrate<- numparticles/numproposals
    hat_p_vec<-sapply(1:N^2,FUN=function(I) sum(ABCdraws_d[I,])/numparticles)

    #if(totnumproposals >= number_sim) return(ABCdraws)
  }
  return(ABCdraws)
}


#'@rdname SMCABCnJRNMM_r1
#'@title SMCABCnJRNMM_r1
#'@description SMC-ABC for time-varying thresholds computed as percentiles of ACCEPTED distances with a fixed budget (number of simulations)
#'@param data observed data
#'@param Pr_cont Matrix of lower/upper bounds for the uniform priors for the continuous parameters
#'@param extra_model List with Number of points to simulate N, final time T, grid, time step h, initial position, vector of diagonal entries of the gamma matrix, vector of diagonal entries of the sigma matrix, Theta0
#'@param ABCthreshold Initial value for the threshold for the guided approach
#'@param minaccept.rate Lowest acceptance rate after which we stop the algorithm
## #'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param numparticles number of particles to keep
#'@param attempt number of iteration
#'@param whichprior choose between 'unif','lognormal' and and 'exp'
#'@param subsamplingby every how many simulated points the observation should be taken. The default=1, i.e., no subsampling
#'@export
#'
SMCABCnJRNMM_r1<- function (data, Pr_cont,extra_model, ABCthreshold, minaccept.rate, numparticles,
                                    attempt, whichprior = 'unif',subsamplingby=1)
{
  N<-extra_model[[1]]
  T<-extra_model[[2]]
  grid<-extra_model[[3]]
  h<-extra_model[[4]]
  startv<-extra_model[[5]]
  dGamma<-extra_model[[6]]
  dSigma<-extra_model[[7]]
  Theta0<-extra_model[[8]]
  ndata<-length(grid)

  #% THE OBSERVED SUMMARIES - Code for the FHN
  summaries_extra<-abc_summaries_extra(data,T,h)
  summobs<-abc_summaries_X(data,T,h,summaries_extra)
  summ_weights<-abc_summaries_weights(summobs,summaries_extra,N)
  #####
  nfreepar_c <- N+2# % number of continuous parameters
  nfreepar_d <- N^2# % number of discrete parameters
  nfreepar <- nfreepar_c+nfreepar_d# % total number of parameters to be inferred (-N in fact)

  ABCdraws <- matrix(0,nrow=nfreepar,ncol=numparticles); #% will store accepted parameters (for 1 iteration only, not all)
  nsummaries <- length(summ_weights)#% the number of summary statistics
  distance_accepted <- matrix(0,nrow=1,ncol=numparticles); #% this is only for the olcm proposal
  lengthmodel<-length(data)
  #% initialization: t is the iteration counter
  t <- 1;
  numproposals0<-0
  number_sim<- numparticles/minaccept.rate
  RES<-foreach(success = 1:numparticles,.combine='rbind',.packages=c('SMCABCnJRNMM')) %dopar% {
    distance<-ABCthreshold+1
    numproposals <- 0;
  #  simsumm_all <- c();

    while(distance >= ABCthreshold & numproposals<number_sim){
      numproposals <- numproposals +1;
      parameters <-  problemprior(-1,1,Pr_cont,whichprior);# % propose from the prior
      Theta0[,1]<- parameters[1:N]
      K<-KmatrixgivenLc(N,parameters[N+1],parameters[N+2])
      Rho<-matrix(parameters[-(1:nfreepar_c)],nrow=N,byrow=T)
      simdata <- model(N, grid, h, startv, dGamma,dSigma, Theta0,Rho,K); #// simulate from the model
      simdata<- simdata[,seq(1,ndata,by=subsamplingby)]
      summaries_Y<-abc_summaries_X(simdata,T,h,summaries_extra)
      simsumm <- abc_summaries(N,summobs,summaries_Y,summaries_extra)
      # xc <- (t(simsumm)-t(summobs)); #// compute
      distance <- abc_distance(simsumm,summ_weights,N)#abc_distance(xc,summ_weights,we)
   #   simsumm_all <-  cbind(simsumm_all,simsumm)
    }
    if(numproposals>=number_sim) {list(numproposals);stop('a particle got stucked');}
    return(list(numproposals,distance,parameters))# return(list(numproposals,distance,parameters,simsumm_all))
  }
  eval_time <- toc()

return(RES)}


#'@rdname abcpilot
#'@title abcpilot
#'@description Reference-table ABC for time-varying thresholds computed as percentiles of ACCEPTED distances with a fixed budget (number of simulations)
#'@param data observed data
#'@param Pr_cont Matrix of lower/upper bounds for the uniform priors for the continuous parameters
#'@param extra_model List with Number of points to simulate N, final time T, grid, time step h, initial position, vector of diagonal entries of the gamma matrix, vector of diagonal entries of the sigma matrix, Theta0
#'@param numsim Number of simulation in the ABC pilot
#'@param whichprior Chosen prior. Here only 'unif' is available
#'@param subsamplingby every how many simulated points the observation should be taken. The default=1, i.e., no subsampling
#'@export
#'
abcpilot<- function (data, Pr_cont,extra_model,numsim, whichprior = 'unif',subsamplingby=1)
{
  N<-extra_model[[1]]
  T<-extra_model[[2]]
  grid<-extra_model[[3]]
  h<-extra_model[[4]]
  startv<-extra_model[[5]]
  dGamma<-extra_model[[6]]
  dSigma<-extra_model[[7]]
  Theta0<-extra_model[[8]]
  ndata<-length(grid)

  #% THE OBSERVED SUMMARIES - Code for the FHN
  summaries_extra<-abc_summaries_extra(data,T,h)
  summobs<-abc_summaries_X(data,T,h,summaries_extra)
  summ_weights<-abc_summaries_weights(summobs,summaries_extra,N)
  #####
  nfreepar_c <- N+2# % number of continuous parameters
  nfreepar_d <- N^2# % number of discrete parameters
  nfreepar <- nfreepar_c+nfreepar_d# % total number of parameters to be inferred (-N in fact)

  #% initialization: t is the iteration counter
  abcdistancepilot<-foreach(success = 1:numsim,.combine='rbind',.packages=c('SMCABCnJRNMM')) %dopar% {
    parameters <-  problemprior(-1,1,Pr_cont,whichprior);# % propose from the prior
    Theta0[,1]<- parameters[1:N]
    K<-KmatrixgivenLc(N,parameters[N+1],parameters[N+2])
    Rho<-matrix(parameters[-(1:nfreepar_c)],nrow=N,byrow=T)
    simdata <- model(N, grid, h, startv, dGamma,dSigma, Theta0,Rho,K); #// simulate from the model
    simdata<- simdata[,seq(1,ndata,by=subsamplingby)]
    summaries_Y<-abc_summaries_X(simdata,T,h,summaries_extra)
    simsumm <- abc_summaries(N,summobs,summaries_Y,summaries_extra)
    return(list(abc_distance(simsumm,summ_weights,N),parameters))
  }
  return(abcdistancepilot)
}


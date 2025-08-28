# SMCABCnJRNMM

R package for the inference of relevant real-valued model parameters and $\{0,1\}$-valued network parameter of N populations of Jansen-and-Rit Neural Mass Models (JRNMM) proposed in 
[1] S. Ditlevsen, M. Tamborrino, I. Tubikanec. Network Inference via Approximate Bayesian Computation. Illustration on a stochastic multi-population Neural Mass Model. To appear in the Annals of Applied Statistics. Preprint at ArXiv: 2306.15787 https://arxiv.org/abs/2306.15787 

The inference is based on an adapted Sequential Monte Carlo Approximate Bayesian computation (SMC-ABC) algorithm, with a proposal sampler accounting for both continuous and binary parameters (adapted in that sense). While applied to the N populations of JRNMM, the method presented in [1], see Algorithm 1, is more general, being a SMC-ABC scheme for network inference (which we denoted nSMCABC).

The R-package is written and maintained by Massimiliano Tamborrino (firstname dot secondname at warwick.ac.uk).

# What can you find in the package
In this package, we provide the code for the nSMC-ABC algorithm (SMC-ABC for network inference) using a Bernoulli proposal for the binary $\{0,1\}$ parameters and a Gaussian kernel proposal for the real-valued ones. The Gaussian proposal is either the canonical version ('standard') or its optimased version ('olcm').  

The main routine is "SMCABCnJRNMM_acceptrate.R", which performs nSMC-ABC with automatically decreasing tolerance levels and uniform priors, and user-specified: 1) Gaussian sampler ('standard' or 'olcm'); 2) minimum acceptance rate below which to stop the algorithm


# How to install the package

* The simplest way is to install the package is via devtools, using devtools::install_github("massimilianotamborrino/SMCABCnJRNMM")

# Output of "SMCABCnJRNMM_acceptrate.R"
The output files (see below) will be saved in the user-specified folder, and will contain information on iteration-stage (t) and attempt. 

Output files: The input "folder" is the name of the folder where you want your results to be stored. The whole path to access to the folder can be given as input. That folder will contain many files upon completion of a run. Here is a description:

- posterior draws: a typical file name ABCdraws_stageX_attemptY.txt This is a matrix of d x N values containing M parameters draws (particles) for each of the d parameters to infer. The X in 'stageX' is the iteration number the draws have been sampled at. The Y in 'attemptY' is the ID of the inference run (since it is possible to run several inference runs on the same dataset, one after the other). 
- ess_stageX_attemptY.txt: Value of the effective sample size (ESS) at the given iteration and attempt.
- evaltime_stageX_attemptY.txt: Number of seconds required to accept N draws at the given iteration and given attempt.
- numproposals_stageX_attemptY.txt: Number of particles proposed (i.e. both the accepted and rejected), at the given iteration and attempt.
- numproposals0_stageX_attemptY.txt: Number of particles drawn from the proposal sampler which have been immediately rejected as out of the prior support, at the given iteration and attempt.
- numproposalsneg_stageX_attemptY.txt: Number of sampled real-valued parameters drawn from the Gaussian proposal sampler which have been immediately rejected as being negative (while all real-valued parameters are positive), at the given iteration and attempt.
- norm_weights_c_stageX_attemptY.txt: Vector of M normalised importance weights associated to each of the M accepted particles at the given iteration and attempt.
- ABCthreshold_stageX_attemptY.txt: Value of the ABC threshold used at the given iteration and attempt.
- Sigma_stageX_attemptY.txt: Covariance matrix used by the Gaussian sampler (either standard or olcm) at the given iteration and attempt.

 

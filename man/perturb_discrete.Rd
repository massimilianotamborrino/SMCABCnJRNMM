% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proposal_kernels.R
\name{perturb_discrete}
\alias{perturb_discrete}
\title{perturb_discrete}
\usage{
perturb_discrete(N, pvec, stay_prob)
}
\arguments{
\item{N}{number of neural populations}

\item{pvec}{N*(N-1)-dimensional probability vector of the binary parameter}

\item{stay_prob}{Probability of not perturbing the sampled binary value}
}
\value{
1 sampled value from the discrete Perturbation kernel for the discrete parameters
}
\description{
Sampling from a Bernoulli distribution, with a further perturbation on the sampled value according to a q_stay probability
}

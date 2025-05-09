% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_summaries.r
\name{abc_summaries_X}
\alias{abc_summaries_X}
\title{abc_summaries_X}
\usage{
abc_summaries_X(X, T, h, summaries_extra)
}
\arguments{
\item{X}{Nxn matrix of observed dataset/path of N populations of JR-NMMs with n points each}

\item{T}{time horizon for path simulation}

\item{h}{step size for path simulation}

\item{summaries_extra}{vector of quantities needed to compute the summary statistics}
}
\value{
summaries (densities, spectral densities, MSCs, ccfs) of a dataset X
}
\description{
Calculations of the ABC summaries of the reference data X
}

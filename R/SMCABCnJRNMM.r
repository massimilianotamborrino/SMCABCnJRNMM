#' SMCABCnJRNMM
#' @description Inference via SMC-ABC for the a network of N stochastics JR-NMMs
#' @docType package
#' @author Massimiliano Tamborrino
#' @import Rcpp StrangSplittingJRNMM mvnfast foreach
#' @useDynLib SMCABCnJRNMM, .registration = TRUE
## #' @importFrom Rcpp sourceCpp
#' @importFrom pracma tic toc isposdef
#' @importFrom Matrix nearPD
## #' @importFrom StrangSplittingJRNMM JRNMM_Splitting_Cpp
## #' @importFrom StrangSplittingJRNMM fast_JRNMM_Splitting_Cpp
## #' @importFrom StrangSplittingJRNMM exponential_matrix_JRNMM_M
## #' @importFrom StrangSplittingJRNMM covariance_matrix_JRNMM_M
#' @importFrom stats runif density ts ts.union ccf median spectrum quantile
#' @importFrom utils write.table
## #' @importFrom IRISSeismic crossSpectrum
## #' @importFrom mvnfast dmvn rmvn
#' @name SMCABCnJRNMM
NULL

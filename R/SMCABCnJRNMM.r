#' SMCABCnJRNMM
#' @description Inference via SMC-ABC for the a network of N stochastics JR-NMMs
#' @docType package
#' @author Massimiliano Tamborrino
#' @import Rcpp StrangSplittingJRNMM mvnfast foreach
#' @useDynLib SMCABCnJRNMM, .registration = TRUE
#' @importFrom pracma tic toc isposdef
#' @importFrom Matrix nearPD
#' @importFrom stats runif density ts ts.union ccf median spectrum quantile
#' @importFrom utils write.table
## #' @importFrom IRISSeismic crossSpectrum
#' @name SMCABCnJRNMM
NULL

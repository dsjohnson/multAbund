#' @title Use Dirichlet Process Prior to Model Multivariate Animal Abundance and Occurance 
#' 
#' @description Functions for implementing the RJMCMC estimation of multivariate abundence and occurence. 
#' These functions model multivariate abundence with latent group 
#' designations for each species. The associations between species are 
#' functions of common group membership.
#' 
#' A vignette is available that illustrates use of the package with simulated abundance data.
#' The vignette can be viewed with the command:
#'  
#' \verb{vignette("simulation_demo", package="multAbund")}
#' 
#' \tabular{ll}{ 
#' Package: \tab multAbund\cr 
#' Type: \tab Package\cr 
#' Version: \tab 0.02\cr 
#' Date: \tab January 12, 2016\cr 
#' License: \tab CC0 \cr 
#' LazyLoad: \tab yes\cr 
#' }
#' 
#' @name multAbund-package
#' @aliases multAbund-package multAbund
#' @docType package
#' @author Devin S. Johnson
#' 
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#'
###' @importFrom Rcpp evalCpp
#' @useDynLib multAbund

NULL

#' @title Use Dirichlet Process Prior to Model Multivariate Animal Abundance and Occurance 
#' 
#' @description Functions for implementing the RJMCMC estimation of multivariate abundence and occurence. 
#' These functions model multivariate abundence with latent group 
#' designations for each species. The associations between species are 
#' functions of common group membership.
#' 
#' \tabular{ll}{ 
#' Package: \tab multAbund\cr 
#' Type: \tab Package\cr 
#' Version: \tab 0.03.9002\cr 
#' Date: \tab March 14, 2022\cr 
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
#' @importFrom Rcpp evalCpp
#' @importFrom stats binomial dgamma glm integrate model.frame model.matrix optim qgamma terms
#' @importFrom utils browseURL 
#' @useDynLib multAbund

NULL

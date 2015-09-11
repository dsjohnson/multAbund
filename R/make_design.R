#' @title Make necessary data matrices for multivariate abundance and occurence RJMCMC functions
#' @param counts A formula that provides the column which contains the counts. Must be of the form \code{~ var - 1} to avoid having a extra column of 1s.
#' @param delta_model A formula decribing the variables for which delta will be the coefficients.
#' @param X_additional A formula giveng additional variables to be used for all groups. 
#' @param sigma_model A formula describing the sigma vector. Must be of the form \code{~factorVariable-1}.
#' @param data Data.frame containing counts and environmental variables
#' 
#' @return 
#' A list containing the following elements:
#' \item{n}{A vector of counts}
#' \item{H}{The environmental variables from which the clustering is based}
#' \item{X}{The design matrix for the global covariates}
#' \item{D}{Design matrix defining sigma_ij}
#' \item{data}{Data frame containing the original data}
#' 
#' @author Devin S. Johnson
#' 
#' @export

make_data_list = function(counts=~1, delta_model, X_additional, sigma_model=~1, data){
  data = with(data, data[order(species, obs), ])
  H_data = model.frame(paste("~", paste(as.character(delta_model)[2], "+ obs")), data)
  H_data = unique(H_data)
  x_form = as.character(delta_model)[2]
  if(!missing(X_additional) ) x_form = paste(x_form, as.character(X_additional)[2], sep=" + ")
  x_form = as.formula(paste("~", x_form))
  return(
    list(
      n = model.matrix(counts, data),
      H=model.matrix(delta_model, H_data),
      X = model.matrix(x_form, data),
      D=model.matrix(sigma_model, data),
      data = data
    )
  )
}

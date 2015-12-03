#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Poisson
#' observations
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param pred_list A named list created in the same way as \code{data_list}, however, 
#' this list will be used for prediction purposes. For example, a user may want to withold
#' data in a cross-validation check.
#' @param initial_vals A named list of initial values for the parameters (see details). 
#' The functions \code{multAbund::sugs} or \code{multAbund::make_inits} can be used to create this
#' list. 
#' @param phi_beta Scale parameter for the prior covariance of \verb{beta}. The covariance 
#' is defined to be \verb{phi_beta^2*solve(X'X)}, where \code{X} is the design matrix 
#' for the covariates associated with \verb{beta}
#' @param mu_beta Prior mean of \verb{beta}.
#' @param phi_omega Scale parameter for omega prior. The omega prior is a half-t
#' distribution.
#' @param df_omega Degrees of freedom for the omega prior distribution. 
#' If \verb{df_omega > 50} then the prior is set to a half-normal distribution rather
#' than a half-t. If \verb{df_omega = 1}, then a half-Cauchy distribution is used.
#' @param a_alpha Shape parameter for the gamma distribution of the \verb{alpha} parameter.
#' @param b_alpha Rate parameter for the gamma distribution of the \verb{alpha} parameter.
#' @param phi_sigma Scale parameter for the half-t prior distribution for \verb{sigma}.
#' @param df_sigma Degrees of freedom for \verb{sigma} prior. Follows the same rules as the omega prior specification.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details Here are some details.
#' @section Vignette:
#' A demonstration using simulated data is available in vignette form and can be accessed 
#' by typing \verb{vignette("simulation_demo", package="multAbund")}. The vignette demo 
#' uses a realistic number of iterations, so, if the user decides to run the associated \verb{R}
#' code it will take some time. 
#' @author Devin S. Johnson
#' @export
mult_abund_pois = function(data_list, pred_list, initial_vals, phi_beta, 
                           mu_beta, phi_omega, df_omega, a_alpha, b_alpha, 
                           phi_sigma, df_sigma, block, begin_group_update, 
                           burn, iter){
  out = mult_abund_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    initial_vals=initial_vals, 
    phi_beta=phi_beta, 
    mu_beta=mu_beta, 
    phi_omega=phi_omega, 
    df_omega=df_omega, 
    a_alpha=a_alpha, 
    b_alpha=b_alpha, 
    phi_sigma=phi_sigma, 
    df_sigma=df_sigma, 
    block=block, 
    begin_group_update=begin_group_update, 
    burn=burn, 
    iter=iter
  )
  return(out)
}


#' @title Perform RJMCMC for posterior sampling of multivariate cluster occurence 
#' model with Bernoulli observations
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param pred_list A named list created in the same way as \code{data_list}, however, 
#' this list will be used for prediction purposes. For example, a user may want to withold
#' data in a cross-validation check.
#' @param initial_vals A named list of initial values for the parameters (see details). 
#' The functions \code{multAbund::sugs} or \code{multAbund::make_inits} can be used to create this
#' list. 
#' @param phi_beta Scale parameter for the prior covariance of \verb{beta}. The covariance 
#' is defined to be \verb{phi_beta^2*solve(X'X)}, where \code{X} is the design matrix 
#' for the covariates associated with \verb{beta}
#' @param mu_beta Prior mean of \verb{beta}.
#' @param phi_omega Scale parameter for omega prior. The omega prior is a half-t
#' distribution.
#' @param df_omega Degrees of freedom for the omega prior distribution. 
#' If \verb{df_omega > 50} then the prior is set to a half-normal distribution rather
#' than a half-t. If \verb{df_omega = 1}, then a half-Cauchy distribution is used.
#' @param a_alpha Shape parameter for the gamma distribution of the \verb{alpha} parameter.
#' @param b_alpha Rate parameter for the gamma distribution of the \verb{alpha} parameter.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details Here are some details.
#' @section Vignette:
#' A demonstration using simulated data is available in vignette form and can be accessed 
#' by typing \verb{vignette("simulation_demo", package="multAbund")}. The vignette demo 
#' uses a realistic number of iterations, so, if the user decides to run the associated \verb{R}
#' code it will take some time. 
#' @author Devin S. Johnson
#' @export
mult_abund_probit = function(data_list, pred_list, initial_vals, phi_beta, 
                             mu_beta, phi_omega, df_omega, a_alpha, b_alpha, 
                             block, begin_group_update, 
                             burn, iter){
  out = mult_occ_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    initial_vals=initial_vals, 
    phi_beta=phi_beta, 
    mu_beta=mu_beta, 
    phi_omega=phi_omega, 
    df_omega=df_omega, 
    a_alpha=a_alpha, 
    b_alpha=b_alpha, 
    block=block, 
    begin_group_update=begin_group_update, 
    burn=burn, 
    iter=iter
  )
  return(out)
}
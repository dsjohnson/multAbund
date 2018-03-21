#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Zero-Inflated Poisson
#' observations
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of Zero-Inflated Poisson count data using 
#' a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
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
#' @param update_omega Logical. Should omega be updated or not?
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
mult_abund_zip = function(data_list, initial_list, prior_list,
                           block, begin_group_update, update_omega=T,
                           burn, iter){
  out = mult_zip_mcmc(
    data_list=data_list, 
    prior_list=prior_list,
    initial_list=initial_list, 
    block=block, 
    begin_group_update=begin_group_update, 
    update_omega=update_omega,
    burn=burn, 
    iter=iter
  )
  return(out)
}

#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Gaussian 
#' observations
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of Gaussian abundance 
#' (typically resulting from a log transformation) data using a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
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
#' @param update_omega Logical. Should omega be updated or not?
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
mult_abund_norm = function(data_list, initial_list, prior_list,
                          block, begin_group_update, update_omega=T,
                          burn, iter){
  out = mult_norm_mcmc(
    data_list=data_list, 
    prior_list=prior_list,
    initial_list=initial_list, 
    block=block, 
    begin_group_update=begin_group_update, 
    update_omega=update_omega,
    burn=burn, 
    iter=iter
  )
  return(out)
}

#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Poisson
#' observations
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of Poisson count data using 
#' a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
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
mult_abund_pois = function(data_list, initial_list, prior_list,
                           block, begin_group_update,
                           burn, iter){
  out = mult_pois_mcmc(
    data_list=data_list, 
    prior_list=prior_list,
    initial_list=initial_list, 
    block=block, 
    begin_group_update=begin_group_update, 
    burn=burn, 
    iter=iter
  )
  return(out)
}


#' @title Perform RJMCMC for posterior sampling of multivariate cluster occurence 
#' model with Bernoulli observations
#' #' @description Fit a Dirichlet Process random effect model for joint species distribution inference of binary occurence data using 
#' a RJMCMC procedure.
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of binary occurence data using 
#' a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
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
#' @author Devin S. Johnson
#' @export
mult_abund_probit = function(data_list, prior_list, initial_list, 
                             block, begin_group_update, 
                             burn, iter){
  out = mult_occ_mcmc(
    data_list=data_list, 
    initial_list=initial_list, 
    prior_list=prior_list,
    block=block, 
    begin_group_update=begin_group_update, 
    burn=burn, 
    iter=iter
  )
  return(out)
}

probit_reg = function(data_list, prior_list, burn, iter){
  inits = glm(data_list$y ~ data_list$X-1, family=binomial("probit"))$coef
  inits = ifelse(is.na(inits), 0, inits)
  out = probit_reg_mcmc(
    data_list=data_list, 
    prior_list=prior_list,
    initial_list=list(beta=inits),
    burn=burn, 
    iter=iter
  )
  return(out)
}

pois_reg = function(data_list, prior_list, initial_list, block, burn, iter){
  inits = glm(data_list$n ~ data_list$X-1, family="poisson")$coef
  inits = ifelse(is.na(inits), 0, inits)
  out = pois_reg_mcmc(
    data_list=data_list, 
    prior_list=prior_list,
    initial_list=list(beta=inits),
    block=block,
    burn=burn, 
    iter=iter
  )
  return(out)
}


zip_reg = function(data_list, prior_list, initial_list, block, burn, iter){
  out = zip_reg_mcmc(
    prior_list=prior_list,
    initial_list=initial_list,
    block=block,
    burn=burn, 
    iter=iter
  )
  return(out)
}

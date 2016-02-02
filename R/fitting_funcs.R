#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Zero-Inflated Poisson
#' observations
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of Poisson count data using 
#' a RJMCMC procedure.
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
mult_abund_zip = function(data_list, pred_list=NULL, initial_list, prior_list,
                           block, begin_group_update, 
                           burn, iter){
  out = mult_zip_mcmc(
    data_list=data_list, 
    pred_list=pred_list,
    prior_list=prior_list,
    initial_list=initial_list, 
    block=block, 
    begin_group_update=begin_group_update, 
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
mult_abund_pois = function(data_list, pred_list=NULL, initial_list, prior_list,
                           block, begin_group_update, 
                           burn, iter){
  out = mult_abund_mcmc(
    data_list=data_list, 
    pred_list=pred_list,
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
#' @author Devin S. Johnson
#' @export
mult_abund_probit = function(data_list, pred_list, prior_list, initial_list, 
                             block, begin_group_update, 
                             burn, iter){
  out = mult_occ_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    initial_list=initial_list, 
    prior_list=prior_list,
    block=block, 
    begin_group_update=begin_group_update, 
    burn=burn, 
    iter=iter
  )
  return(out)
}

#' @title Probit regression via MCMC
#' @description Fit a probit regression model via MCMC using data augmentation.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param pred_list A named list created in the same way as \code{data_list}, however, 
#' this list will be used for prediction purposes. For example, a user may want to withold
#' data in a cross-validation check.
#' @param prior_list A list of prior parameters with the names and elements: 
#' (1) \verb{Sigma_beta_inv}, the inverse covariance matrix of the normal prior for beta, and 
#' (2) \verb{mu_beta}, the prior mean of beta.
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @author Devin S. Johnson
#' @export
probit_reg = function(data_list, pred_list, burn, iter){
  inits = glm(data_list$y ~ data_list$X-1, family=binomial("probit"))$coef
  inits = ifelse(is.na(inits), 0, inits)
  out = probit_reg_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    prior_list-prior_list,
    initial_list=list(beta=inits),
    burn=burn, 
    iter=iter
  )
  return(out)
}

#' @title Poisson regression via MCMC
#' @description Fit a Poisson regression model with normal randome effect for modeling overdispersion.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param pred_list A named list created in the same way as \code{data_list}, however, 
#' this list will be used for prediction purposes. For example, a user may want to withold
#' data in a cross-validation check.
#' @param prior_list A list of prior parameters with the names and elements: 
#' (1) \verb{Sigma_beta_inv}, the inverse covariance matrix of the normal prior for beta,
#' (2) \verb{mu_beta}, the prior mean of beta, 
#' (3) \verb{phi_sigma} scale parameter for the half-t prior distribution for sigma,
#' (4) \verb{df_sigma} the degrees of freedom for sigma prior. Follows the same rules as the omega prior specification.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @author Devin S. Johnson
#' @export
pois_reg = function(data_list, pred_list, prior_list, initial_list, block, burn, iter){
  inits = glm(data_list$n ~ data_list$X-1, family="poisson")$coef
  inits = ifelse(is.na(inits), 0, inits)
  out = pois_reg_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    prior_list=prior_list,
    initial_list=list(beta=inits),
    block=block,
    burn=burn, 
    iter=iter
  )
  return(out)
}

#' @title ZIP regression via MCMC
#' @description Fit a Zero-Inflated Poisson (ZIP) regression model with normal randome effect for modeling overdispersion.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param pred_list A named list created in the same way as \code{data_list}, however, 
#' this list will be used for prediction purposes. For example, a user may want to withold
#' data in a cross-validation check.
#' @param prior_list A list of prior parameters with the names and elements: 
#' (1) \verb{Sigma_beta_inv}, the inverse covariance matrix of the normal prior for beta,
#' (2) \verb{mu_beta}, the prior mean of beta, 
#' (3) \verb{phi_sigma} scale parameter for the half-t prior distribution for sigma,
#' (4) \verb{df_sigma} the degrees of freedom for sigma prior. Follows the same rules as the omega prior specification.
#' @param initial_list named list of initial values. 
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @author Devin S. Johnson
#' @importFrom pscl zeroinfl
#' @export
zip_reg = function(data_list, pred_list, prior_list, initial_list, block, burn, iter){
  out = zip_reg_mcmc(
    data_list=data_list, 
    pred_list=pred_list, 
    prior_list=prior_list,
    initial_list=initial_list,
    block=block,
    burn=burn, 
    iter=iter
  )
  return(out)
}

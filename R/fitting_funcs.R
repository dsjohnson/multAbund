#' @title Perform RJMCMC for posterior sampling of multivariate cluster abundance model with Zero-Inflated Poisson
#' observations
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of Zero-Inflated Poisson count data using 
#' a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
#' The functions \code{multAbund::sugs} or \code{multAbund::make_inits} can be used to create this
#' list. 
#' @param prior_list A named list of prior distribution parameters. See details.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param update_omega Logical. Should omega be updated or not?
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details The \code{prior_list} argument needs to contain the following items:
#' \itemize{
#' \item a_alpha numeric. shape parameter for gamma prior on alpha,
#' \item b_alpha numeric. scale parameter for gamma prior on alpha,
#' \item Sigma_beta_inv numeric matrix. precision matrix for global regression coefficients, beta
#' \item mu_beta numeric vector. prior mean of beta parameters
#' \item phi_omega numeric. scale for half t/normal prior on delta variance parameter (omega),
#' \item df_omega numeric. degrees of freedom for half-t prior on omega (df_omega>=50 means a half-normal will be assumed)
#' \item phi_sigma numeric. scale parameter for half-t/normal prior for sigma parameters
#' \item df_sigma numeric. degrees of freedom for sigma prior (>=50 implies half-normal will be used)
#' }
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
#' (typically resulting from a log transformation) data using a RJMCMC procedure. THIS IS CURRENTLY A
#' TEST FUNCTION THAT IS NOT FULLY DEVELOPED YET! IT USES A DIFFERENT VARIANCE MODEL FOR DELTA PARAMETERS.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
#' The functions \code{multAbund::sugs} or \code{multAbund::make_inits} can be used to create this
#' list. 
#' @param prior_list A named list of prior distribution parameters. See details.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param update_omega Logical. Should omega be updated or not?
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details The \code{prior_list} argument needs to contain the following items:
#' \itemize{
#' \item a_alpha numeric. shape parameter for gamma prior on alpha,
#' \item b_alpha numeric. scale parameter for gamma prior on alpha,
#' \item Sigma_beta_inv numeric matrix. precision matrix for global regression coefficients, beta
#' \item mu_beta numeric vector. prior mean of beta parameters
#' \item phi_omega numeric. scale for half t/normal prior on delta variance parameter (omega),
#' \item df_omega numeric. degrees of freedom for half-t prior on omega (df_omega>=50 means a half-normal will be assumed)
#' \item phi_sigma numeric. scale parameter for half-t/normal prior for sigma parameters
#' \item df_sigma numeric. degrees of freedom for sigma prior (>=50 implies half-normal will be used)
#' }
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
  om_var = attr(terms(data_list$delta_model), "term.labels")
  if(attr(terms(data_list$delta_model), "intercept") == 1) om_var = c("(Intercept)", om_var) 
  colnames(out$omega) = om_var
  colnames(out$beta) = colnames(data_list$X)
  colnames(out$delta_bar) = do.call(paste, expand.grid(paste('/',colnames(data_list$H)), levels(factor(data_list$data$species)))[,2:1])
  colnames(out$sigma) = colnames(data_list$D)
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
#' @param prior_list A named list of prior parameters. See details.
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details The \code{prior_list} argument needs to contain the following items:
#' \itemize{
#' \item a_alpha numeric. shape parameter for gamma prior on alpha,
#' \item b_alpha numeric. scale parameter for gamma prior on alpha,
#' \item Sigma_beta_inv numeric matrix. precision matrix for global regression coefficients, beta
#' \item mu_beta numeric vector. prior mean of beta parameters
#' \item phi_omega numeric. scale for half t/normal prior on delta variance parameter (omega),
#' \item df_omega numeric. degrees of freedom for half-t prior on omega (df_omega>=50 means a half-normal will be assumed)
#' \item phi_sigma numeric. scale parameter for half-t/normal prior for sigma parameters
#' \item df_sigma numeric. degrees of freedom for sigma prior (>=50 implies half-normal will be used)
#' }
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
#' @description Fit a Dirichlet Process random effect model for joint species distribution inference of binary occurence data using 
#' a RJMCMC procedure.
#' @param data_list A named list created of data items created from user data with the 
#' function \code{multAbund::make_data_list}. This data will be used for model fitting.
#' @param initial_list A named list of initial values for the parameters (see details). 
#' The functions \code{multAbund::sugs} or \code{multAbund::make_inits} can be used to create this
#' list.
#' @param prior_list A named list of prior parameters. See details.  
#' @param block Number of iterations between Metropolis proposal adaptation.
#' @param begin_group_update The iteration at with the group clusters begin updating.
#' The RJMCMC often performs better when the chain is allowed to sample only the parameters 
#' before the groups begin updating. 
#' @param update_omega logical. Should omega be updated during the MCMC or remain fixed.
#' @param burn Number of burnin iterations that are discarded.
#' @param iter Number of iterations retained for posterior inference.
#' @details The \code{prior_list} argument needs to contain the following items:
#' \itemize{
#' \item a_alpha numeric. shape parameter for gamma prior on alpha,
#' \item b_alpha numeric. scale parameter for gamma prior on alpha,
#' \item Sigma_beta_inv numeric matrix. precision matrix for global regression coefficients, beta
#' \item mu_beta numeric vector. prior mean of beta parameters
#' \item phi_omega numeric. scale for half t/normal prior on delta variance parameter (omega),
#' \item df_omega numeric. degrees of freedom for half-t prior on omega (df_omega>=50 means a half-normal will be assumed)
#' \item phi_sigma numeric. scale parameter for half-t/normal prior for sigma parameters
#' \item df_sigma numeric. degrees of freedom for sigma prior (>=50 implies half-normal will be used)
#' }
#' @author Devin S. Johnson
#' @export
mult_abund_probit = function(data_list, prior_list, initial_list, 
                             block, begin_group_update, update_omega,
                             burn, iter){
  out = mult_occ_mcmc(
    data_list=data_list, 
    initial_list=initial_list, 
    prior_list=prior_list,
    block=block, 
    begin_group_update=begin_group_update, 
    update_omega=update_omega,
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

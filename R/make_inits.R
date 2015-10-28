find_alpha=function(kappa, n){
  foo = function(alpha, kappa, n){
    return((kappa-alpha*sum(1/(alpha+c(1:n)-1)))^2)
  }
  return(optimize(foo, c(0,100), kappa=kappa, n=n)$minimum)
}

#' @title Compute intial values for the abundance / occurence RJMCMC using K-means clustering
#' @param data_list A list created by the make_data_list function.
#' @param log_alpha An intial guess at the log of the Dirichlet process parameter. Defaults to 0.
#' @param phi_beta The prior variance parameter for the Normal prior for the fixed effect beta parameters
#' @param mu_beta The prior mean for beta
#' @param phi_omega The scale parameter for the half-t prior distribution for omega
#' @param df_omega Degrees-of-freedom parmeter for omega half-t prior distribution. 
#' @return A list with the following elements:
#' \item{beta}{A vector with the initial beta value}
#' \item{delta}{A vector with the initial delta values}
#' \item{groups}{Initial group assignments for each species}
#' \item{omega}{Inital value for omega}
#' \item{sigma}{Vector of inital values for sigma. Currently this is a vector filled with exp(-10). i.e., very little overdispersion.}
#' \item{log_alpha}{Initial value of log alpha for the Dirichlet process prior}
#' @author Devin S. Johnson
#' @importFrom mvtnorm dmvnorm
#' @importFrom NbClust NbClust
#' @export

make_inits = function(
  data_list,
  log_alpha=0,
  phi_beta=10, 
  mu_beta, 
  phi_omega = 1, 
  df_omega = 1
){
  data = data_list$data
  data$species = factor(data$species)
  num_spec = length(levels(data$species))
  if(!is.null(data_list$n)){
    resp = data_list$n
    family = "poisson"
  } else{
    data$resp = data_list$y
    family = binomial(link = "probit")
  }
  H = data_list$H
  X = data_list$X
  groups = 1:num_spec
  K_pi = kronecker(diag(num_spec), H)
  # fit1 = glm(n ~ K_pi - 1, family="poisson")
  #beta = solve(crossprod(X), crossprod(X, log(n+1)))
  fit = glm(resp ~ X + K_pi - 1, family=family)
  delta = fit$coef[-c(1:ncol(X))]
  beta = fit$coef[1:ncol(X)]
  delta[is.na(delta)] = 0
  delta_mat = matrix(delta, ncol=ncol(H), byrow=TRUE)
  pdf(file = NULL)
  sink("/dev/null")
  nc = suppressWarnings(NbClust(delta_mat, min.nc = 2, max.nc=num_spec-2, method="kmeans"))
  sink()
  dev.off()
  num_groups = max(nc$Best.partition)
  groups_init = kmeans(delta_mat, num_groups)
  groups = groups_init$cluster
  delta = as.vector(t(groups_init$centers))
  #delta_mat = sweep(delta_mat, 2, apply(delta_mat, 2, mean), "-")
  omega=optimize(
    f=function(x, H, df_omega){
      -sum(dmvnorm(delta_mat, rep(0,ncol(H)), x^2*solve(crossprod(H)),log=TRUE)) - 
        dt(x/phi_omega, df_omega, log=TRUE) - log(phi_omega)
    }, 
    interval=c(0,1000), H=H, df_omega=df_omega
  )$minimum
  sigma = rep(1.0e-4, ncol(data_list$D))
  alpha = find_alpha(max(groups), num_spec)
  
  return(
    list(
      beta=beta,
      delta=delta,
      groups=groups,
      omega=omega,
      sigma=sigma, 
      log_alpha=log(alpha)
    ))
  
}

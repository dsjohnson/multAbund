find_alpha=function(kappa, n){
  foo = function(alpha, kappa, n){
    return((kappa-alpha*sum(1/(alpha+c(1:n)-1)))^2)
  }
  return(optimize(foo, c(0,100), kappa=kappa, n=n)$minimum)
}

#' @title Compute intial values for the abundance / occurence RJMCMC using K-means clustering
#' @param data_list A list created by the make_data_list function.
#' @param phi_omega The scale parameter for the half-t prior distribution for omega
#' @param df_omega Degrees-of-freedom parmeter for omega half-t prior distribution. 
#' @param num_groups An inital guess at the number of groups
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

make_inits = function(
  data_list,
  phi_omega = 1, 
  df_omega = 1,
  num_groups
){
  data = data_list$data
  data$species = factor(data$species)
  num_spec = length(levels(data$species))
  if(!is.null(data_list$n)){
    resp = data_list$n
    family = "poisson"
  } else{
    resp = data_list$y
    family = binomial(link = "probit")
  }
  H = data_list$H
  X = data_list$X
  groups = 1:num_spec
  K_pi = kronecker(diag(num_spec), H)
  fit = glm(resp ~ X + K_pi - 1, family=family)
  delta = fit$coef[-c(1:ncol(X))]
  beta = fit$coef[1:ncol(X)]
  delta[is.na(delta)] = 0
  delta_mat = matrix(delta, ncol=ncol(H), byrow=TRUE)
  omega=optimize(
    f=function(x, H, df_omega){
      -sum(dmvnorm(delta_mat, rep(0,ncol(H)), x^2*solve(crossprod(H)),log=TRUE)) - 
        dt(x/phi_omega, df_omega, log=TRUE) - log(phi_omega)
    }, 
    interval=c(0,1000), H=H, df_omega=df_omega
  )$minimum
  if(missing(num_groups)){
    pdf(file = NULL)
    sink("/dev/null")
    nc = suppressWarnings(NbClust(delta_mat, min.nc = 2, max.nc=num_spec-2, method="kmeans"))
    sink()
    dev.off()
    num_groups = max(nc$Best.partition)
  }
  if(num_groups==num_spec){
    delta=as.vector(t(delta_mat))
    groups=c(1:num_spec)
  } else if(num_groups==1){
    delta = colMeans(delta_mat)
    groups=rep(1,num_spec)
  } else{
    groups_init = kmeans(delta_mat, num_groups)
    groups = groups_init$cluster
    delta = as.vector(t(groups_init$centers))
  }
  #delta_mat = sweep(delta_mat, 2, apply(delta_mat, 2, mean), "-")
  omega=optimize(
    f=function(x, H, df_omega){
      -sum(dmvnorm(delta_mat, rep(0,ncol(H)), x^2*solve(crossprod(H)),log=TRUE)) - 
        dt(x/phi_omega, df_omega, log=TRUE) - log(phi_omega)
    }, 
    interval=c(0,1000), H=H, df_omega=df_omega
  )$minimum
  # omega=2.85
  sigma = rep(1.0e-4, ncol(data_list$D))
  alpha = find_alpha(max(groups), num_spec)
  
  return(
    list(
      beta=beta,
      delta=delta,
      groups=as.integer(groups),
      omega=omega,
      sigma=sigma, 
      log_alpha=log(alpha)
    ))
  
}

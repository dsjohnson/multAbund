find_alpha=function(kappa, n){
  foo = function(alpha, kappa, n){
    return((kappa-alpha*sum(1/(alpha+c(1:n)-1)))^2)
  }
  return(optimize(foo, c(0,100), kappa=kappa, n=n)$minimum)
}

#' @title Compute intial values for the abundance / occurence RJMCMC using an approximation of the SUGS algorithm
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
#' @references Wang, L. and Dunson, D.B. (2011) Fast Bayesian Inference in Dirichlet Process Mixture Models. Journal of Computational and Graphical Statistics, 20:196-216.
#' @importFrom mvtnorm dmvnorm
#' @export

sugs = function(
  data_list,
  prior_list,
  log_alpha
#   log_alpha=0,
#   phi_beta=10, 
#   mu_beta, 
#   phi_omega = 1, 
#   df_omega = 1
){
  data = data_list$data
  if(!is.null(data_list$n)){
    resp = data_list$n
    family = "poisson"
  } else{
    resp = data_list$y
    family = binomial(link = "probit")
  }
  phi_omega=prior_list$phi_omega
  df_omega=prior_list$df_omega
  Sigma_beta_inv = prior_list$Sigma_beta_inv
  Sigma_beta = solve(Sigma_beta_inv)
  mu_beta=prior_list$mu_beta
  H = data_list$H
  X = data_list$X
  groups = rep(1,max(data_list$data$species))
  shuffle = sample(1:max(data_list$data$species), max(data_list$data$species))
  K_pi = kronecker(diag(max(data$species)), H)
  fit = glm(resp ~ X + K_pi - 1, family=family)
  delta = fit$coef[-c(1:ncol(X))]
  delta[is.na(delta)]=0
  beta = fit$coef[1:ncol(X)]
  delta_mat = matrix(delta, ncol=ncol(H), byrow=TRUE)
  delta_mat = sweep(delta_mat, 2, apply(delta_mat, 2, mean), "-")
  omega=optimize(
    f=function(x, H, df_omega){
      -sum(dmvnorm(delta_mat, rep(0,ncol(H)), x^2*solve(crossprod(H)),log=TRUE)) - 
        dt(x/phi_omega, df_omega, log=TRUE) - log(phi_omega)
    }, 
    interval=c(0,1000), H=H, df_omega=df_omega
  )$minimum
  for(i in 2:max(data$species)){
    resp_tmp = resp[data$species%in%shuffle[1:i]]
    X_tmp = X[data$species%in%shuffle[1:i],]
    idx = (apply(X_tmp, 2, var)!=0) | colMeans(X_tmp)==1
    X_tmp = X_tmp[,idx]
    Sigma_beta_tmp = Sigma_beta[idx,idx]
    ln_I = rep(0,max(groups)+1)
    ln_table=c(log(as.integer(table(groups[shuffle[1:(i-1)]]))), log_alpha)
    for(h in 1:(max(groups)+1)){
      groups[shuffle[i]] = h
      if(all(groups[shuffle[1:i]]==1)){
        d = ncol(X_tmp)
        fit = glm(resp_tmp ~ X_tmp-1, family=family)
        ln_I[h] = (d/2) * log(2*pi) + 
          0.5*log(det(vcov(fit))) + 
          as.double(logLik(fit)) +
          dmvnorm(coefficients(fit), mu_beta[idx], Sigma_beta_tmp, log=TRUE) + 
          # dmvnorm(rep(0,ncol(H)), rep(0,ncol(H)), omega^2 * solve(crossprod(H)), log=TRUE) +
          ln_table[h]
      } else{
        C_pi = model.matrix(~factor(groups[shuffle[1:i]])-1)
        K_pi = kronecker(C_pi, H)
        d = ncol(X_tmp) + ncol(K_pi)
        fit = glm(resp_tmp ~ X_tmp + K_pi - 1, family=family)
        delta = fit$coef[-c(1:ncol(X_tmp))]
        beta = fit$coef[1:ncol(X_tmp)]
        delta[is.na(delta)] = 0
        delta_mat = matrix(delta, ncol=ncol(H), byrow=TRUE)
        delta_mat = sweep(delta_mat, 2, apply(delta_mat, 2, mean), "-")
        ln_I[h] = (d/2) * log(2*pi) + 
          0.5*log(det(vcov(fit))) +
          as.double(logLik(fit)) + 
          dmvnorm(beta, mu_beta[idx], Sigma_beta_tmp, log=TRUE) + 
          sum(dmvnorm(delta_mat, rep(0,ncol(H)), omega^2 * solve(crossprod(H)), log=TRUE)) +
          ln_table[h]
      }
    } # end table loop
    groups[shuffle[i]] = which(ln_I == max(ln_I))
  } # end species loop
  C_pi = model.matrix(~factor(groups)-1)
  K_pi = kronecker(C_pi, H)
  d = ncol(X_tmp) + ncol(K_pi)
  fit = glm(n ~ X + K_pi - 1, family="poisson")
  sigma = rep(1.0e-4, ncol(data_list$D))
  delta = fit$coef[-c(1:ncol(X))]
  beta = fit$coef[1:ncol(X)]
  delta[is.na(delta)] = 0
  delta_mat = matrix(delta, ncol=ncol(H), byrow=TRUE)
  delta_mat = sweep(delta_mat, 2, apply(delta_mat, 2, mean), "-")
  omega=optimize(
    f=function(x, H, df_omega){
      -sum(dmvnorm(delta_mat, rep(0,ncol(H)), x^2*solve(crossprod(H)),log=TRUE)) - 
        dt(x, df_omega, log=TRUE)
    }, 
    interval=c(0,100), H=H, df_omega=df_omega
  )$minimum
  alpha=find_alpha(max(groups), nrow(C_pi))
  ln_I = (d/2) * log(2*pi) + 
    0.5*log(det(vcov(fit))) +
    as.double(logLik(fit)) + 
    dmvnorm(beta, mu_beta[idx], Sigma_beta_tmp, log=TRUE) + 
    sum(dmvnorm(delta_mat, rep(0,ncol(H)), omega^2 * solve(crossprod(H)), log=TRUE))
    
  
  
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

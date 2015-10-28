
foo = function(x, n, k, mu, phi){
  exp(k*log(x) + lgamma(x) - lgamma(x+n) + dgamma(x, mu, phi, log=TRUE))  #
}

#' @importFrom copula Stirling1
group_dist = function(n, mu=0, phi=1.5){
  out=rep(NA, n)
  for(k in 1:n){
    out[k] = abs(copula::Stirling1(n,k))*integrate(foo, 0, Inf, n=n, k=k, phi=phi, mu=mu)$value
  }
  return(out/sum(out))
}

#' @title Obtain parameters for gamma prior such that the distribution of the number of groups is equal to the \code{target} argument
#' @description This function uses the method of Dorazio (2009) to calculate the parameters of a gamma prior on the Dirichlet process such that the induced distribution function
#' of the number of clusters is equal to 1/k.
#' @param n An integer representing the maximum number of groups. Typically this is the number of individuals. 
#' @param target The desired prior probability mass function for the number of groups. This function must return a 
#' vector of length \code{n} which contains the prior probabilities for each group size.  
#' @author Devin S. Johnson
#' @references Dorazio, R. M. (2009) On selecting a prior for the precision parameter of Dirichlet process mixture models. Journal of Statistical Planning and Inference 139:3384-3390.
#' @export

get_opt_alpha_prior = function(n, target=NULL){
  if(is.null(target)) target = function(n){x = 1/c(1:n); return(x/sum(x))}
  ab = optim(c(0.3,0.02), function(par, n, target){sum(abs(group_dist(n,par[1],par[2]) - target(n))^2)}, 
             n=n, target=target, lower=c(0.1,0.0005), method="L-BFGS-B")$par
  return(
    list(
      ab=ab,
      pmf = group_dist(n, ab[1], ab[2])
    )
  )
}


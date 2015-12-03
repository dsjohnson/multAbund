
foo = function(x, n, k, a, b){
  abs(copula::Stirling1(n,k)) * exp(k*log(x) + lgamma(x) - lgamma(x+n) + dgamma(x, a, b, log=TRUE))  #
}

#' @importFrom copula Stirling1
group_dist = function(n, a, b){
  out=rep(NA, n)
  ul = qgamma(0.999, a, b)
  ll = qgamma(0.001, a, b)
  for(k in 1:n){
    out[k] = integrate(foo, 0, Inf, n=n, k=k, a=a, b=b, stop.on.error = FALSE)$value
  }
  return(out/sum(out))
}

dkl = function(g, tgt){
  sum(tgt*(log(tgt)-log(g)))
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
  ab = exp(optim(c(log(0.4),log(0.01)), 
             function(par, n, target){
               return(dkl(group_dist(n,exp(par[1]),exp(par[2])), target(n)))
             }, n=n, target=target)$par)
  return(
    list(
      ab=ab,
      pmf = group_dist(n, ab[1], ab[2])
    )
  )
}


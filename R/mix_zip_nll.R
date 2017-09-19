logSumExp = function(x){
  xmax <- which.max(x)
  log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
}

ln_zip = function(n, gamma, lambda){
  log(gamma*(n==0) + (1-gamma)*dpois(n,lambda))
}


mix_zip_nll = function(par, n, species, X, H, M, n_groups){
  if(is.null(X)) nX = 0 else nX = ncol(X)
  npar = c(nX, ncol(H)*n_groups, ncol(M), n_groups-1)
  par_list = split(par, rep(c('beta','delta','logit_gamma','logit_pi'), npar))
  delta = matrix(par_list$delta, ncol=n_groups)
  Xb = X%*%par_list$beta
  Hdelta = H %*% delta
  gamma = plogis(M%*%par_list$logit_gamma)
  pi = exp(c(0, par_list$logit_pi))/sum(exp(c(0, par_list$logit_pi)))
  ll = matrix(NA,length(n), n_groups)
  for(k in 1:n_groups){
    ll[,k] = ln_zip(n, gamma, exp(Xb + Hdelta[,k]))
  }
  ll = aggregate(ll, list(species), sum)[,-1]
  ll = t(apply(ll, 1, function(x,ln_pi) x+ln_pi, ln_pi = log(pi)))
  return(-sum(apply(ll, 1, logSumExp)))
}

mle_mult_abund_pois = function(counts, species, beta_form=NULL, 
                               delta_form=~1, gamma_form=~1, n_groups, ln_prior_func=NULL,
                               data, init=NULL, nrep=1, ...){
  if(!is.null(beta_form)) X = model.matrix(beta_form, data) else X=NULL
  H = model.matrix(delta_form, data)
  M = model.matrix(gamma_form, data)
  n = data[,counts]
  species = data[,species]
  if(is.null(ln_prior_func)) ln_prior = function(par){return(0)}
  if (is.null(init)) {
    beta = solve(crossprod(X),crossprod(X, log(n+0.5)))
    delta = solve(crossprod(H),crossprod(H, log(n+0.5)-X%*%beta))
    pi = c(n_groups:1)/sum(c(n_groups:1))
    pi = log(pi/pi[1])[-1]
    par = c(
      beta,
      rep(delta, n_groups),
      rep(0, ncol(M)),
      pi
    )
    obj_func = function(par, n, species, X, H, M, n_groups){
      return(mix_zip_nll(par, n, species, X, H, M, n_groups) - ln_prior_func(par))
    }
    fit_list = vector(mode = 'list', nrep)
    for(rr in 1:nrep){
      init = optim(par, obj_func, n=n, species=species, X=X, H=H, M=M, n_groups=n_groups,
                  method='SANN', control = list(maxit=2000))
      message(rr, ')', ' initial -lnl = ', init$value)
      fit_list[[rr]] = optim(init$par, obj_func, n=n, species=species, X=X, H=H, M=M, n_groups=n_groups,...)
      if(fit_list[[rr]]$convergence==0){
        message('     converged -lnl = ', fit_list[[rr]]$value)
      } else{
        message('     ** not converged **')
      }
    }
    best = which.min(sapply(fit_list, function(x) x$value))
    out = fit_list[[best]]
  }
  else {
    par = init
    out = optim(par, mix_zip_nll, n=n, species=species, X=X, H=H, M=M, n_groups=n_groups,...)
  }
  return(out)
}
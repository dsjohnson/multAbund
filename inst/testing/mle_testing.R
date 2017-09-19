delta_form = ~V1+V2+V3
beta_form = ~V4-1
gamma_form=~1
n_groups=5
counts = "count"
species = "species"
data = cnt_dat

delta = sweep(matrix(delta_pi, 5, 4, byrow=T), 2, beta[-5], "+")

init = c(0.5, as.vector(t(delta)), -10, -0.3364722, -0.5596158, -0.8472979, -1.9459101)

alpha = 6
Sig = matrix(trigamma(alpha/n_groups), n_groups-1, n_groups-1)
diag(Sig) = rep(2*trigamma(alpha/n_groups))
dir_prior = function(par){
  return(dmvnorm(tail(par,4), rep(0,4), Sig, log = T))
}

fit = mle_mult_abund_pois(counts, species, beta_form, delta_form, 
                          gamma_form, n_groups, ln_prior_func=NULL, 
                          data=data, nrep=3, #init=init,
                          control=list(maxit=20000))


mix_zip_nll(par, n, species, X, H, M, n_groups)

debug(mle_mult_abund_pois)
debug(mix_zip_nll)

undebug(mle_mult_abund_pois)
undebug(mix_zip_nll)

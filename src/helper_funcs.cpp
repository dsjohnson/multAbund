// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include <progress.hpp>

using namespace Rcpp;
using namespace arma;

int sample_du(arma::vec ppp){
  arma::vec cdf = cumsum(ppp/sum(ppp));
  double U = Rcpp::as<double>(Rcpp::runif(1));
  int out = 1;
  if(U<= cdf[0]) return(out);
  else
  {
    for(int i=1; i<ppp.n_elem; i++){ 
      if(U <= cdf[i]){
        out = i+1;
        return(out);
      }
    }
    return(out);
  }
}

arma::vec armaNorm(int n){
  NumericVector x = rnorm(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}

arma::vec armaU(int n){
  NumericVector x = runif(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}

arma::vec GCN(const arma::mat& M, const arma::vec& m){
  arma::vec out = solve(M,m) + solve(chol(M), armaNorm(m.n_elem));
  return out;
}

arma::vec ln_Pois(const arma::vec& k, const arma::vec& ln_lambda){
  return k%ln_lambda - exp(ln_lambda);
}

double ln_norm(const arma::vec& x, const arma::vec& sigma_sq){
  NumericVector out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out[i] = R::dnorm(x(i), 0, sqrt(sigma_sq(i)),1);
  }
  return sum(out);
}

double ln_norm_2(const arma::vec& x, const double& sigma_sq){
  NumericVector out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out[i] = R::dnorm(x(i), 0, sqrt(sigma_sq),1);
  }
  return sum(out);
}

double ln_mvnorm(const arma::vec& x, const arma::mat& Sigma_inv){
  return -0.5*x.n_elem*log(2*PI) + 0.5*log(det(Sigma_inv)) - 0.5*as_scalar(x.t()*Sigma_inv*x);
}

double ln_cauchy(const arma::vec& x, const double& scale){
  arma::vec out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out(i) = R::dcauchy(x(i), 0, scale, 1);
  }
  return sum(out);
}

double ln_t(const arma::vec& x, const double& scale, const double& df){
  arma::vec out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    if(df<50){
      out(i) = R::dt(x(i)/scale, df, 1) - log(scale);
    } else {
      out(i) = R::dnorm(x(i), 0, scale, 1);
    }
  }
  return sum(out);
}

double ln_logis(const arma::vec& x, const double& mu, const double& scale){
  arma::vec out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out(i) = R::dlogis(x(i), mu, scale, 1);
  }
  return sum(out);
}

double ln_logistic_beta(const arma::vec& x, const double& a, const double& b){
  arma::vec out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out(i) = R::dbeta(R::plogis(x(i),0,1,1,0), a, b, 1) + x(i) - 2*log(1+exp(x(i)));
  }
  return sum(out);
}


double ln_t_2(const double& x, const double& scale, const double& df){
  double out;
  if(df<50){
    out = R::dt(x/scale, df, 1) - log(scale);
  } else {
    out = R::dnorm(x, 0, scale, 1);
  }
  return out;
}

double ln_zip(const int& x, const double& p, const double& lb){
  if(x==0){
    return(log(p + (1-p)*R::dpois(x, lb, 0)));
  } else{
    return(log(1-p) + R::dpois(x, lb, 1));
  }
}

double ln_zip_vec(const arma::vec& x, const arma::vec& p, const arma::vec& lb){
  NumericVector out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    if(x(i)==0){
      out[i] = log(p(i) + (1-p(i))*R::dpois(x(i), lb(i), 0));
    } else{
      out[i] = log(1-p(i)) + R::dpois(x(i), lb(i), 1);
    }
  }
  return(sum(out));
}

arma::vec logit_inv(const arma::vec& x){
  arma::vec out(x.n_elem);
  for(int i=0; i<out.n_elem; i++){
    out(i) = R::plogis(x(i),0,1,1,0);
  }
  return(out);
}

arma::vec logit(const arma::vec& p){
  arma::vec out(p.n_elem);
  for(int i=0; i<out.n_elem; i++){
    out(i) = R::qlogis(p(i),0,1,1,0);
  }
  return(out);
}




double ln_crp(const double& log_alpha, const arma::mat& C_pi){
  arma::vec groups = sum(C_pi).t();
  double alpha = exp(log_alpha);
  double out =  lgamma(alpha) - lgamma(alpha + C_pi.n_rows) + C_pi.n_cols*log_alpha;
  for(int i=0; i<groups.n_elem; i++) out += lgamma(groups(i));
  return out;
}

arma::mat LtoC(const arma::mat& L){
  arma::vec G(L.n_rows, fill::ones);
  int kappa;
  for(int i=1; i<L.n_rows; i++){
    kappa = max(G);
    if(L(i,i)==1){
      G(i) = kappa+1;
    } else{
      for(int j=0; j<i; j++){
        if(L(i,j)==1) G(i) = G(j);
      }
    }
  }
  kappa = max(G);
  arma::mat C(L.n_rows, kappa, fill::zeros);
  for(int i=0; i<L.n_rows; i++){
    C(i,G(i)-1) = 1;
  }
  return C;
}

arma::mat rmult(const arma::vec& sigma, const arma::mat& X){
  arma::mat C = X;
  for(int i=0; i<C.n_cols; i++){
    C.col(i) %= sigma;
  }
  return C;
}

arma::vec arma_pnorm(const arma::vec& x){
  arma::vec out(x.n_elem);
  for(int i=0; i<x.n_elem; i++){
    out(i) = R::pnorm(x(i),0,1,1,0);
  }
  return out;
}

arma::vec arma_rtruncnorm(const arma::vec& mean, const arma::vec& a, const arma::vec& b){
  arma::vec out(mean.n_elem);
  double truncLo;
  double truncUp;
  for(int i=0; i<mean.n_elem; i++){
    truncLo = R::pnorm(a(i), mean(i), 1, 1, 0);
    truncUp = R::pnorm(b(i), mean(i), 1, 1, 0);
    out(i) = R::qnorm(R::runif(truncLo, truncUp), mean(i), 1, 1, 0);
  }
  return out;
}

arma::vec arma_rbern(const arma::vec& p){
  arma::vec out(p.n_elem);
  for(int i=0; i<p.n_elem; i++){
    out(i) = R::rbinom(1,p(i));
  }
  return(out);
}

arma::vec arma_rpois(const arma::vec& lam){
  arma::vec out(lam.n_elem);
  for(int i=0; i<lam.n_elem; i++){
    out(i) = R::rpois(lam(i));
  }
  return(out);
}

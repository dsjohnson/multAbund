// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include <progress.hpp>

using namespace Rcpp;
using namespace arma;

// Function defs..
int sample_du(arma::vec ppp);
arma::vec armaNorm(int n);
arma::vec armaU(int n);
arma::vec GCN(const arma::mat& V, const arma::vec& v);
double ln_norm(const arma::vec& x, const arma::vec& sigma_sq);
double ln_norm_2(const arma::vec& x, const double& sigma_sq);
double ln_mvnorm(const arma::vec& x, const arma::mat& Sigma_inv);
double ln_cauchy(const arma::vec& x, const double& scale);
double ln_t(const arma::vec& x, const double& scale, const double& df);
double ln_t_2(const double& x, const double& scale, const double& df);
double ln_crp(const double& log_alpha, const arma::mat& C_pi);
arma::mat LtoC(const arma::mat& L);
arma::mat rmult(const arma::vec& sigma, const arma::mat& X);
arma::vec arma_pnorm(const arma::vec& x);
arma::vec arma_rtruncnorm(const arma::vec& mean, const arma::vec& a, const arma::vec& b);
arma::vec arma_rbern(const arma::vec& p);


// [[Rcpp::export]]
List probit_reg_mcmc(
    const Rcpp::List& data_list,
    const Rcpp::List& prior_list,
    const Rcpp::List& initial_list,
    const int& burn, 
    const int& iter
) {
  
  //matrices for fitting
  arma::vec y = as<arma::vec>(data_list["y"]);
  arma::mat X = as<arma::mat>(data_list["X"]);
//   arma::mat X_pred;
//   if(!Rf_isNull(pred_list)){
//     X_pred = as<arma::mat>(pred_list["X"]);
//   } else{
//     X_pred = X;
//   }
  arma::uvec obs_idx = find_finite(y);
  X = X.rows(obs_idx);
  y = y.elem(obs_idx);
  const int I = X.n_rows;
  const int p = X.n_cols;

  // #beta items
  arma::vec beta = as<arma::vec>(initial_list["beta"]);;
  arma::mat Sigma_beta_inv = as<arma::mat>(prior_list["Sigma_beta_inv"]);
  arma::vec mu_beta = as<arma::vec>(prior_list["mu_beta"]);
  arma::mat  V_beta_inv(p,p);
  arma::vec v_beta(p);
  arma::mat beta_store(iter, p);
  
  // #z items
  arma::vec z_prop(I);
  arma::mat z_store(iter, I);
  arma::vec r_z(I);
  arma::vec z(y.n_elem);
  arma::vec mu_z = X*beta;
  arma::vec a(y.n_elem);
  arma::vec b(y.n_elem);
  for(int j=0;j<y.n_elem;j++){
    if(y(j)==0){
      a(j) = -datum::inf;
      b(j) = 0;
    } else if(y(j)==1){
      a(j) = 0;
      b(j) = datum::inf;
    } else{
      a(j) = -datum::inf;
      b(j) = datum::inf;
    }
  }
  
  // Rcout << "H" << endl;
  
  // other quantities 
  // arma::mat pred_store(iter, X_pred.n_rows);
  
  //Begin MCMC
  Progress prog(iter+burn, true);
  
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z = X*beta;
    z = arma_rtruncnorm(mu_z,a,b);
    if(i>=burn) z_store.row(i-burn) = z.t();
    
    // update beta
    V_beta_inv = X.t()*X + Sigma_beta_inv;
    v_beta = X.t()*z + Sigma_beta_inv*mu_beta;
    beta = GCN(V_beta_inv, v_beta);
    if(i>=burn) beta_store.row(i-burn) = beta.t();
    
    // make prediction
//     if(i>=burn){
//       pred_store.row(i-burn) = arma_rbern(arma_pnorm(X_pred*beta)).t();
//     }
    
    prog.increment();
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("z") = z_store,
    Rcpp::Named("beta") = beta_store//,
    // Rcpp::Named("pred")=pred_store
  );  
}
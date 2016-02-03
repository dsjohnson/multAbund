
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include <progress.hpp>

using namespace Rcpp;
using namespace arma;


arma::vec armaNorm(int n);
arma::vec armaU(int n);
arma::vec GCN(const arma::mat& V, const arma::vec& v);
double ln_norm(const arma::vec& x, const arma::vec& sigma_sq);
double ln_norm_2(const arma::vec& x, const double& sigma_sq);
double ln_mvnorm(const arma::vec& x, const arma::mat& Sigma_inv);
double ln_t(const arma::vec& x, const double& scale, const double& df);
double ln_t_2(const double& x, const double& scale, const double& df);
arma::mat rmult(const arma::vec& sigma, const arma::mat& X);
double ln_zip(const int& x, const double& p, const double& lb);
double ln_zip_vec(const arma::vec& x, const arma::vec& p, const arma::vec& lb);
arma::vec logit_inv(const arma::vec& x);
arma::vec logit(const arma::vec& p);
arma::vec arma_rbern(const arma::vec& p);
double ln_logis(const arma::vec& x, const double& mu, const double& scale);
arma::vec arma_rpois(const arma::vec& lam);

// [[Rcpp::export]]
List zip_reg_mcmc(
    const Rcpp::List& data_list,
    const Rcpp::List& pred_list,
    const Rcpp::List& prior_list,
    const Rcpp::List& initial_list,
    const int& block, 
    const int& burn, 
    const int& iter
) {
  
  //matrices for fitting
  arma::vec n = as<arma::vec>(data_list["n"]);
  arma::mat X = as<arma::mat>(data_list["X"]);
  arma::mat D = as<arma::mat>(data_list["D"]);
  arma::mat M = as<arma::mat>(data_list["M"]);
  arma::mat X_pred;
  arma::mat D_pred;
  arma::mat M_pred;
  X_pred = as<arma::mat>(pred_list["X"]);
  D_pred = as<arma::mat>(pred_list["D"]);
  M_pred = as<arma::mat>(pred_list["M"]);
  arma::uvec obs_idx = find_finite(n);
  X = X.rows(obs_idx);
  n = n.elem(obs_idx);
  D = D.rows(obs_idx);
  M = M.rows(obs_idx);
  const int I = X.n_rows;
  const int p = X.n_cols;
  
  // Rcout << "Prelim OK" << endl;
  
  
  // #sigma items
  double df_sigma = as<double>(prior_list["df_sigma"]);
  double phi_sigma = as<double>(prior_list["phi_sigma"]);
  arma::vec log_sigma = -10*ones<vec>(D.n_cols);
  arma::vec log_sigma_prop = log_sigma;
  double MHR_sigma;
  arma::mat log_sigma_store(iter+burn, log_sigma.n_elem);
  arma::vec jump_sigma(iter+burn, fill::zeros);
  double r_sigma;
  double tune_log_sigma = 2.4*2.4/log_sigma.n_elem;
  arma::mat pv_log_sigma = 0.1*eye(log_sigma.n_elem,log_sigma.n_elem);
  arma::mat L_log_sigma = chol(pv_log_sigma).t();
  arma::vec sigma2_z = exp(2*D*log_sigma) + 1.0E-8;
  arma::vec sigma2_z_prop = sigma2_z;
  
  // Rcout << "sigma prelim OK" << endl;
  
  // #beta items
  arma::mat Sigma_beta_inv = as<arma::mat>(prior_list["Sigma_beta_inv"]);
  arma::vec mu_beta = as<arma::vec>(prior_list["mu_beta"]);
  arma::mat  V_beta_inv(p,p);
  arma::vec v_beta(p);
  arma::vec beta(X.n_cols);
  beta = as<arma::vec>(initial_list("beta")); 
  arma::mat beta_store(iter, p);
  
  // gamma items 
  double phi_lg = as<double>(prior_list["phi_logit_gamma"]);
  double df_lg = as<double>(prior_list["df_logit_gamma"]);
  arma::vec lg = as<arma::vec>(initial_list("logit_gamma"));
  arma::vec lg_prop = lg;
  double MHR_lg;
  arma::mat lg_store(iter+burn, lg.n_elem);
  arma::vec jump_lg(iter+burn, fill::zeros);
  double r_lg;
  double tune_lg = 2.4*2.4/lg.n_elem;
  arma::mat pv_lg = 0.1*eye(lg.n_elem,lg.n_elem);
  arma::mat L_lg = chol(pv_lg).t();
  arma::vec gamma = logit_inv(M*lg);
  arma::vec gamma_prop = gamma;
  arma::vec gamma_pred = logit_inv(M_pred*lg);
  
  
  
  // Rcout << "beta prelim OK" << endl;
  
  // #z items
  arma::vec z_prop(I);
  arma::vec U(I);
  arma::vec MHR_z(I);
  arma::mat jump_z(iter+burn, I, fill::zeros);
  arma::uvec jump;
  arma::vec tune_z(I);
  tune_z.fill(2.4*2.4);
  arma::vec pv_z(I, fill::ones);
  arma::mat z_store(iter+burn, I);
  arma::vec r_z(I);
  arma::uvec i_uvec(1);
  arma::mat tune_store(burn+iter, I);
  arma::vec z = log(n+0.5);
  arma::vec mu_z = X*beta; 
  
  // Rcout << "z prelim ok" << endl;
  
  // other quantities 
  arma::mat pred_store(iter, X_pred.n_rows);
  arma::vec mu_z_pred(X_pred.n_rows);
  arma::vec z_pred(X.n_rows);
  arma::vec sigma2_z_pred(D_pred.n_rows);
  arma::vec zrs(I, fill::zeros);
  arma::vec n_pred(n.n_elem, fill::zeros);
  
  
  // Rcout << "Storage ok" << endl;
  
  //Begin MCMC
  Progress prog(iter+burn, true);
  
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z = X*beta;
    z_prop = z + sqrt(tune_z)%sqrt(pv_z)%armaNorm(I);
    for(int j=0; j<I; j++){
      z_prop(j) = z(j) + sqrt(tune_z(j))*sqrt(pv_z(j))*R::norm_rand();
      MHR_z(j) = exp(
        ln_zip(n(j), gamma(j), exp(z_prop(j))) + R::dnorm(z_prop(j), mu_z(j), sqrt(sigma2_z(j)), 1)
        - ln_zip(n(j), gamma(j), exp(z(j))) - R::dnorm(z(j),mu_z(j),sqrt(sigma2_z(j)),1)
      );
    }
    jump = find(armaU(I)<=MHR_z);
    z.elem(jump) = z_prop.elem(jump);
    i_uvec(0)=i;
    jump_z(i_uvec,jump) = ones<rowvec>(jump.n_elem);
    z_store.row(i) = z.t();
    
    // adapt z MH tuning parameter
    if(i>0 & i%block==0){
      r_z = mean(jump_z.submat(i-block, 0, i, I-1)).t();
      tune_z = exp(log(tune_z) + pow(i/block,-0.5)*(r_z-0.234));
      pv_z = pv_z + pow(i/block,-0.5)*(var(z_store.submat(i-block, 0, i, I-1)).t() - pv_z);
    }
    tune_store.row(i) = (sqrt(tune_z)%sqrt(pv_z)).t();
    
    // Rcout << "z updated" << endl;
    
    // update beta
    V_beta_inv = X.t()*rmult(1/sigma2_z,X) + Sigma_beta_inv;
    v_beta = X.t()*(z/sigma2_z) + Sigma_beta_inv*mu_beta;
    beta = GCN(V_beta_inv, v_beta);
    mu_z = X*beta;
    if(i>=burn) beta_store.row(i-burn) = beta.t();
    
    // Rcout << "beta updated" << endl;
    
    //update logit gamma
    lg_prop = lg + sqrt(tune_lg)*L_lg*armaNorm(lg.n_elem);
    gamma_prop = logit_inv(M*lg_prop);
    MHR_sigma = exp(
      ln_zip_vec(n, gamma_prop, exp(z)) + sum(ln_logis(lg_prop, 0,1)) //sum(ln_t(lg_prop, phi_lg, df_lg))
      - ln_zip_vec(n, gamma, exp(z)) - sum(ln_logis(lg, 0,1)) //sum(ln_t(lg, phi_lg, df_lg))
    );
    if(R::runif(0,1) <= MHR_sigma){
      lg = lg_prop;
      gamma = gamma_prop;
      jump_lg(i) = 1;
    }
    lg_store.row(i) = lg.t();
    
    // adapt logit(gamma) MH tuning parameter
    if(i>0 & i%block==0){
      r_lg = mean(jump_lg.subvec(i-block, i));
      tune_lg = exp(log(tune_lg) + pow(i/block,-0.5)*(r_lg-0.234));
      pv_lg = pv_lg + sqrt(block/i)*(cov(lg_store.rows(i-block, i)) - pv_lg);
    }
    
    // update sigma
    log_sigma_prop = log_sigma + sqrt(tune_log_sigma)*L_log_sigma*armaNorm(log_sigma.n_elem);
    sigma2_z_prop = exp(2*D*log_sigma_prop) + 1.0E-8;
    MHR_sigma = exp(
      ln_norm(z-mu_z, sigma2_z_prop) + ln_t(exp(log_sigma_prop), phi_sigma, df_sigma) + sum(log_sigma_prop)
      - ln_norm(z-mu_z, sigma2_z) - ln_t(exp(log_sigma), phi_sigma, df_sigma) - sum(log_sigma)
    );
    if(R::runif(0,1) <= MHR_sigma){
      log_sigma = log_sigma_prop;
      sigma2_z =sigma2_z_prop;
      jump_sigma(i) = 1;
    }
    log_sigma_store.row(i) = log_sigma.t();
    
    // adapt log(sigma) MH tuning parameter
    if(i>0 & i%block==0){
      r_sigma = mean(jump_sigma.subvec(i-block, i));
      tune_log_sigma = exp(log(tune_log_sigma) + pow(i/block,-0.5)*(r_sigma-0.234));
      pv_log_sigma = pv_log_sigma + sqrt(block/i)*(cov(log_sigma_store.rows(i-block, i)) - pv_log_sigma);
    }
    
    // Rcout << "sigma updated" << endl;
    
    // make prediction
    if(i>=burn){
      mu_z_pred = X_pred*beta; 
      sigma2_z_pred = exp(D_pred*log_sigma);
      z_pred = mu_z_pred + sigma2_z_pred%armaNorm(X_pred.n_rows);
      gamma_pred = logit_inv(M_pred*lg);
      n_pred = arma_rpois(exp(z_pred));
      pred_store.row(i-burn) = (arma_rbern(gamma_pred)%n_pred).t();
    }
    
    // Rcout << "prediction made" << endl;
    
    prog.increment();
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("z") = z_store.rows(burn, burn+iter-1),
    Rcpp::Named("beta") = beta_store,
    Rcpp::Named("logit_gamma") = lg_store,
    Rcpp::Named("sigma")=exp(log_sigma_store.rows(burn,iter+burn-1)),
    Rcpp::Named("pred")=pred_store
  );  
}
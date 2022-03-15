// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include <progress.hpp>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;
//using namespace R;

// Function defs..
int sample_du(arma::vec ppp);
arma::vec armaNorm(int n);
arma::vec armaU(int n);
arma::vec GCN(const arma::mat& V, const arma::vec& v);
arma::vec ln_Pois(const arma::vec& k, const arma::vec& ln_lambda);
double ln_norm(const arma::vec& x, const arma::vec& sigma_sq);
double ln_norm_2(const arma::vec& x, const double& sigma_sq);
double ln_mvnorm(const arma::vec& x, const arma::mat& Sigma_inv);
double ln_cauchy(const arma::vec& x, const double& scale);
double ln_t(const arma::vec& x, const double& scale, const double& df);
double ln_t_2(const double& x, const double& scale, const double& df);
double ln_crp(const double& log_alpha, const arma::mat& C_pi);
arma::mat LtoC(const arma::mat& L);
arma::mat rmult(const arma::vec& sigma, const arma::mat& X);
arma::vec arma_rpois(const arma::vec& lam);

// [[Rcpp::export]]
List mult_norm_mcmc(
    const Rcpp::List& data_list,
    const Rcpp::List& initial_list, 
    const Rcpp::List& prior_list, 
    const int& block, 
    const int& begin_group_update,
    const bool& update_omega,
    const int& burn, 
    const int& iter
) {
  
  //matrices for fitting
  arma::vec z_orig = as<arma::vec>(data_list["z"]);
  arma::vec z = z_orig;
  arma::mat H = as<arma::mat>(data_list["H"]);
  arma::mat X = as<arma::mat>(data_list["X"]);
  arma::mat D = as<arma::mat>(data_list["D"]);
  arma::mat G = as<arma::mat>(data_list["G"]);
  
  const int J = H.n_rows;
  const int I = X.n_rows/H.n_rows;
  const int q = H.n_cols;
  const int p = X.n_cols;
  const arma::mat I_J(J,J, fill::eye);
  const arma::mat I_I(I,I, fill::eye);
  
  // Rcout << "Matrices loaded" << endl;
  
  // pi items 
  IntegerVector initial_groups = as<IntegerVector>(initial_list["groups"]);
  int kappa_pi = max(initial_groups);
  arma::mat C_pi(I, kappa_pi, fill::zeros); 
  for(int i=0; i<I; i++){C_pi(i,initial_groups(i)-1)=1;}
  arma::mat Prox = C_pi*C_pi.t();
  arma::cube prox_store(I,I,iter);
  // arma::cube L_store(I,I,iter);
  arma::mat L(I,I,fill::zeros);
  L(0,0) = 1;
  for(int i=1; i<I; i++){
    L(i,min(find(Prox.row(i)==1))) = 1;
  }
  arma::vec ln_link_probs;
  arma::rowvec zeros_row(I, fill::zeros);
  arma::mat K_pi = kron(C_pi, H);
  arma::vec res(I*J);
  arma::vec delta_pi_hat;
  arma::vec kappa_pi_store(iter);
  
  // Rcout << "Pi items defined" << endl;
  
  // #sigma items
  double phi_sigma = as<double>(prior_list["phi_sigma"]);
  double df_sigma = as<double>(prior_list["df_sigma"]);
  arma::vec log_sigma = log(as<arma::vec>(initial_list["sigma"]));
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
  
  // Rcout << "C" << endl;
  
  // #omega items
  double phi_omega = as<double>(prior_list["phi_omega"]);
  double df_omega = as<double>(prior_list["df_omega"]);
  arma::vec log_omega_scale = log(as<arma::vec>(prior_list["omega_scale"]));
  arma::mat HtH = H.t()*H;
  arma::mat HtH_inv = inv_sympd(HtH);
  arma::vec log_omega = log(as<arma::vec>(initial_list["omega"]));
  arma::vec log_omega_prop = log_omega;
  arma::mat Omega_inv = diagmat(1/(exp(2*G*log_omega +2*log_omega_scale) + 1.0e-8));
  arma::mat Omega_inv_prop = Omega_inv;
  arma::mat log_omega_store(iter+burn, log_omega.n_elem);
  double MHR_omega;
  arma::vec jump_omega(iter+burn, fill::zeros);
  double r_omega;
  double tune_log_omega = 2.4*2.4/log_omega.n_elem;
  arma::mat pv_log_omega = 0.1*eye(log_omega.n_elem,log_omega.n_elem);
  arma::mat L_log_omega = chol(pv_log_omega).t();
  
  // Rcout << "D" << endl;
  
  // #alpha items
  double a_alpha = as<double>(prior_list["a_alpha"]);
  double b_alpha = as<double>(prior_list["b_alpha"]);
  double log_alpha = as<double>(initial_list["log_alpha"]);
  double log_alpha_prop = 0;
  arma::vec log_alpha_store(iter+burn);
  double MHR_alpha;
  arma::vec jump_alpha(iter+burn, fill::zeros);
  double r_alpha;
  double tune_log_alpha = 2.4*2.4;
  double pv_log_alpha = 1;
  
  // Rcout << "E" << endl;
  
  // #beta items
  arma::mat Sigma_beta_inv = as<arma::mat>(prior_list["Sigma_beta_inv"]);
  arma::vec mu_beta = as<arma::vec>(prior_list["mu_beta"]);
  arma::mat  V_beta_inv(p,p);
  arma::vec v_beta(p);
  arma::vec beta = as<arma::vec>(initial_list["beta"]);
  arma::mat beta_store(iter, p);
  
  // Rcout << "F" << endl;
  
  // #delta_pi items
  arma::mat Sigma_delta_pi = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv.i());
  arma::mat Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv);
  arma::mat Sigma_delta_pi_inv_prop;
  arma::mat V_delta_pi_inv(kappa_pi*q,kappa_pi*q);
  arma::mat V_delta_pi(kappa_pi*q,kappa_pi*q);
  arma::vec v_delta_pi(kappa_pi*q);
  arma::vec delta_pi = as<arma::vec>(initial_list["delta"]);
  arma::mat delta_bar_store(iter, I*q);
  arma::mat A = kron(ones<mat>(kappa_pi,1), eye<mat>(q,q)).t();
  arma::mat to_bar = kron(C_pi, eye<mat>(q,q));
  
  // Rcout << "G" << endl;
  
  // #z items
  arma::vec mu_z = X*beta + K_pi*delta_pi;
  arma::mat z_store(iter+burn, I*J);
  
  //Begin MCMC
  Progress prog(iter+burn, true);
  
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z = X*beta + K_pi*delta_pi;
    for(int j=0; j<I*J; j++){
      if(!is_finite(z_orig(j))){
        z(j) = R::rnorm(mu_z(j),sqrt(sigma2_z(j)));
      }
    }
    z_store.row(i) = z.t();
    
    // Rcout << "z updated" << endl;
    
    // update beta
    V_beta_inv = X.t()*rmult(1/sigma2_z,X) + Sigma_beta_inv;
    v_beta = X.t()*((z - K_pi*delta_pi)/sigma2_z) + Sigma_beta_inv*mu_beta;
    beta = GCN(V_beta_inv, v_beta);
    if(i>=burn) beta_store.row(i-burn) = beta.t();
    
    // Rcout << "beta updated" << endl;
    
    if(i > begin_group_update){
      // update C_pi matrix
      res = z-X*beta;
      for(int k=1; k<I; k++){
        ln_link_probs.set_size(k+1);
        ln_link_probs.zeros();
        for(int link=0; link<=k; link++){
          L.row(k) = zeros_row;
          L(k,link) = 1;
          C_pi = LtoC(L);
          kappa_pi = C_pi.n_cols;
          K_pi = kron(C_pi, H);
          Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv);
          V_delta_pi_inv = K_pi.t()*rmult(1/sigma2_z,K_pi) + Sigma_delta_pi_inv; 
          delta_pi_hat = solve(V_delta_pi_inv, K_pi.t()*(res/sigma2_z));
          ln_link_probs(link) = (kappa_pi*q/2)*log(2*M_PI) 
            - 0.5*log(det(V_delta_pi_inv))
            + ln_norm(res-K_pi*delta_pi_hat, sigma2_z) 
            + ln_mvnorm(delta_pi_hat, Sigma_delta_pi_inv);
        }
        ln_link_probs -= max(ln_link_probs);
        ln_link_probs(k) += log_alpha;
        L.row(k) = zeros_row;
        L(k,sample_du(exp(ln_link_probs))-1) = 1;
      }
      C_pi = LtoC(L);
      kappa_pi = C_pi.n_cols;
      K_pi = kron(C_pi, H);
      to_bar = kron(C_pi, eye<mat>(q,q));
      A = kron(ones<mat>(kappa_pi,1), eye<mat>(q,q)).t();
    }
    if(i>=burn){
      prox_store.slice(i-burn)=C_pi*C_pi.t();
      // L_store.slice(i-burn) = L;
      kappa_pi_store(i-burn) = kappa_pi;
    }
    // Rcout << "C_pi updated" << endl;
    
    // update delta_pi
    Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv);
    V_delta_pi_inv = K_pi.t()*rmult(1/sigma2_z, K_pi) + Sigma_delta_pi_inv; 
    v_delta_pi = K_pi.t()*((z-X*beta)/sigma2_z);
    delta_pi = GCN(V_delta_pi_inv, v_delta_pi);
    V_delta_pi = inv(V_delta_pi_inv);
    delta_pi = delta_pi - V_delta_pi*A.t()*inv(A*V_delta_pi*A.t())*A*delta_pi; 
    if(i>=burn){
      delta_bar_store.row(i-burn) = (to_bar*delta_pi).t(); 
    }
    
    // Rcout << "delta updated" << endl;
    
    // update omega
    if(update_omega){
      log_omega_prop = log_omega + sqrt(tune_log_omega)*L_log_omega*armaNorm(log_omega.n_elem);
      Omega_inv_prop = diagmat(1/(exp(2*G*log_omega_prop + 2*log_omega_scale) + 1.0e-8));
      Sigma_delta_pi_inv_prop = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv_prop);
      Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), Omega_inv);
     if(kappa_pi>1){
        MHR_omega = exp(
          ln_mvnorm(delta_pi, Sigma_delta_pi_inv_prop) 
        + ln_t(exp(log_omega_prop), phi_omega, df_omega) + sum(log_omega_prop)
        - ln_mvnorm(delta_pi, Sigma_delta_pi_inv) 
        - ln_t(exp(log_omega), phi_omega, df_omega) - sum(log_omega)
        );
      } else{
        MHR_omega = exp(
          ln_t(exp(log_omega_prop), phi_omega, df_omega) + sum(log_omega_prop)
        - ln_t(exp(log_omega), phi_omega, df_omega) - sum(log_omega)
        );
      }
      
      // Rcout << MHR_omega << endl;
      
      if(R::runif(0,1) <= MHR_omega){
        log_omega = log_omega_prop;
        Omega_inv = Omega_inv_prop;
        jump_omega(i) = 1;
      }
      log_omega_store.row(i) = log_omega.t();
      
      // adapt log(omega) MH tuning parameter
      if(i>block & i%block==0 & i<= begin_group_update){
        r_omega = mean(jump_omega.subvec(i-block, i));
        tune_log_omega = exp(log(tune_log_omega) + 2*pow(2, 0.25)*pow(i/block,-0.25)*(r_omega-0.234));
        pv_log_omega = pv_log_omega + pow(i/block,-0.25)*(cov(log_omega_store.rows(i-block, i)) - pv_log_omega);
        // Rcout << r_omega<< endl << endl << pv_log_omega << endl;
        L_log_omega = chol(pv_log_omega).t();
        // Rcout << "r_omega = " << r_omega << "   tune_log_omega = " << tune_log_omega << endl;
      }
      
      // Rcout << "omega updated" << endl;
    }
    // update alpha
    log_alpha_prop = log_alpha + sqrt(tune_log_alpha*pv_log_alpha)*R::rnorm(0,1);
    MHR_alpha = exp(
      ln_crp(log_alpha_prop, C_pi) + R::dgamma(exp(log_alpha_prop), a_alpha, 1/b_alpha, 1) + log_alpha_prop
      -ln_crp(log_alpha, C_pi) - R::dgamma(exp(log_alpha), a_alpha, 1/b_alpha, 1) - log_alpha
    );
    if(R::runif(0,1) <= MHR_alpha){
      log_alpha = log_alpha_prop;
      jump_alpha(i) = 1;
    }
    log_alpha_store(i) = log_alpha;
    
    // adapt log(alpha) MH tuning parameter
    if(i>block & i%block==0){
      r_alpha = mean(jump_alpha.subvec(i-block, i));
      tune_log_alpha = exp(log(tune_log_alpha) + 2*pow(2, 0.25)*pow(i/block,-0.25)*(r_alpha-0.234));
      pv_log_alpha = pv_log_alpha + pow(i/block,-0.25)*(var(log_alpha_store.subvec(i-block, i)) - pv_log_alpha);
    }
    
    // Rcout << "alpha updated" << endl;
    
    // update sigma
    mu_z = X*beta + K_pi*delta_pi;
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
    if(i>block & i%block==0){
      r_sigma = mean(jump_sigma.subvec(i-block, i));
      tune_log_sigma = exp(log(tune_log_sigma) + 2*pow(2, 0.25)*pow(i/block,-0.25)*(r_sigma-0.234));
      pv_log_sigma = pv_log_sigma + pow(i/block,-0.25)*(cov(log_sigma_store.rows(i-block, i)) - pv_log_sigma);
      L_log_sigma = chol(pv_log_sigma).t();
    }
    
    // Rcout << "sigma updated" << endl;
    
    prog.increment();
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("z") = z_store.rows(burn, burn+iter-1),
    Rcpp::Named("beta") = beta_store,
    Rcpp::Named("delta_bar") = delta_bar_store, 
    Rcpp::Named("prox") = prox_store,
    Rcpp::Named("kappa_pi") = kappa_pi_store,
    Rcpp::Named("omega")=exp(log_omega_store.rows(burn,iter+burn-1)),
    Rcpp::Named("alpha")=exp(log_alpha_store(span(burn, burn+iter-1))),
    Rcpp::Named("sigma")=exp(log_sigma_store.rows(burn,iter+burn-1))//,
    // Rcpp::Named("pred")=pred_store
  );  
}
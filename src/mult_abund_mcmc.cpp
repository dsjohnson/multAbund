// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "RcppArmadillo.h"
#include <progress.hpp>

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
List mult_pois_mcmc(
    const Rcpp::List& data_list,
    const Rcpp::List& pred_list,
    const Rcpp::List& initial_list, 
    const Rcpp::List& prior_list, 
    const int& block, 
    const int& begin_group_update,
    const int& burn, 
    const int& iter
) {
  
  //matrices for fitting
  arma::vec n = as<arma::vec>(data_list["n"]);
  arma::mat H = as<arma::mat>(data_list["H"]);
  arma::mat X = as<arma::mat>(data_list["X"]);
  arma::mat D = as<arma::mat>(data_list["D"]);
  arma::mat H_pred;
  arma::mat X_pred;
  arma::mat D_pred;
  if(!Rf_isNull(pred_list)){
    H_pred = as<arma::mat>(pred_list["H"]);
    X_pred = as<arma::mat>(pred_list["X"]);
    D_pred = as<arma::mat>(pred_list["D"]);
  } else{
    H_pred = H;
    X_pred = X;
    D_pred = D;
  }
  
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
  arma::mat HtH = H.t()*H;
  arma::mat HtH_inv = inv_sympd(HtH);
  double log_omega = log(as<double>(initial_list["omega"]));
  double log_omega_prop = 0;
  arma::vec log_omega_store(iter+burn);
  double MHR_omega;
  arma::vec jump_omega(iter+burn, fill::zeros);
  double r_omega;
  double tune_log_omega = 2.4*2.4;
  double pv_log_omega = 1;
  
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
  arma::mat Sigma_delta_pi = kron(eye<mat>(kappa_pi,kappa_pi), exp(2*log_omega)*HtH_inv);
  arma::mat Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), exp(-2*log_omega)*HtH);
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
  arma::vec z_prop(I*J);
  arma::vec U(I*J);
  arma::vec MHR_z(I*J);
  arma::mat jump_z(iter+burn, I*J, fill::zeros);
  arma::uvec jump;
  arma::vec tune_z(I*J);
  tune_z.fill(2.4*2.4);
  arma::vec pv_z(I*J, fill::ones);
  arma::mat z_store(iter+burn, I*J);
  arma::vec r_z(I*J);
  arma::uvec i_uvec(1);
  arma::mat tune_store(burn+iter, I*J);
  arma::vec z = log(n+1);
  arma::vec mu_z = X*beta + K_pi*delta_pi;
  
  // Rcout << "H" << endl;
  
  // other quantities 
  arma::mat pred_store(iter, X_pred.n_rows);
  arma::mat K_pi_pred;
  arma::vec mu_z_pred(X_pred.n_rows);
  arma::vec z_pred(X.n_rows);
  arma::vec sigma2_z_pred(D_pred.n_rows);
  
  // Rcout << "I" << endl;
  
  //Begin MCMC
  Progress prog(iter+burn, true);
  
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z = X*beta + K_pi*delta_pi;
    z_prop = z + sqrt(tune_z)%sqrt(pv_z)%armaNorm(I*J);
    for(int j=0; j<I*J; j++){
      if(is_finite(n(j))){
        z_prop(j) = z(j) + sqrt(tune_z(j))*sqrt(pv_z(j))*R::norm_rand();
        MHR_z(j) = exp(
          R::dpois(n(j), exp(z_prop(j)), 1) + R::dnorm(z_prop(j), mu_z(j), sqrt(sigma2_z(j)), 1)
          - R::dpois(n(j), exp(z(j)), 1) - R::dnorm(z(j),mu_z(j),sqrt(sigma2_z(j)),1)
        );
      } else{
        z_prop(j) = R::rnorm(mu_z(j),sqrt(sigma2_z(j)));
        MHR_z(j) = 1.0;
      }
    }
    jump = find(armaU(I*J)<=MHR_z);
    z.elem(jump) = z_prop.elem(jump);
    i_uvec(0)=i;
    jump_z(i_uvec,jump) = ones<rowvec>(jump.n_elem);
    z_store.row(i) = z.t();
    
    // Rcout << "z updated" << endl;
    
    // adapt z MH tuning parameter
    if(i>0 & i%block==0){
      r_z = mean(jump_z.submat(i-block, 0, i, I*J-1)).t();
      tune_z = exp(log(tune_z) + pow(i/block,-0.5)*(r_z-0.234));
      pv_z = pv_z + pow(i/block,-0.5)*(var(z_store.submat(i-block, 0, i, I*J-1)).t() - pv_z);
    }
    tune_store.row(i) = (sqrt(tune_z)%sqrt(pv_z)).t();
    
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
          Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), exp(-2*log_omega)*HtH);
          V_delta_pi_inv = K_pi.t()*rmult(1/sigma2_z,K_pi) + Sigma_delta_pi_inv; 
          delta_pi_hat = solve(V_delta_pi_inv, K_pi.t()*(res/sigma2_z));
          ln_link_probs(link) = (kappa_pi*q/2)*log(2*PI) 
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
    Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), exp(-2*log_omega)*HtH);
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
    log_omega_prop = log_omega + sqrt(tune_log_omega*pv_log_omega)*R::rnorm(0,1);
    Sigma_delta_pi_inv_prop = kron(eye<mat>(kappa_pi,kappa_pi), exp(-2*log_omega_prop)*HtH);
    Sigma_delta_pi_inv = kron(eye<mat>(kappa_pi,kappa_pi), exp(-2*log_omega)*HtH);
    if(kappa_pi>1){
      MHR_omega = exp(
        ln_mvnorm(delta_pi, Sigma_delta_pi_inv_prop) 
      + ln_t_2(exp(log_omega_prop), phi_omega, df_omega) + log_omega_prop
      - ln_mvnorm(delta_pi, Sigma_delta_pi_inv) 
      - ln_t_2(exp(log_omega), phi_omega, df_omega) - log_omega
      );
    } else{
      MHR_omega = exp(
        ln_t_2(exp(log_omega_prop), phi_omega, df_omega) + log_omega_prop
      - ln_t_2(exp(log_omega), phi_omega, df_omega) - log_omega
      );
    }
    if(R::runif(0,1) <= MHR_omega){
      log_omega = log_omega_prop;
      jump_omega(i) = 1;
    }
    log_omega_store.row(i) = log_omega;
  // Rcout << "omega updated" << endl;
    
    // adapt log(omega) MH tuning parameter
    if(i>0 & i%block==0 & i<= begin_group_update){
      r_omega = mean(jump_omega.subvec(i-block, i));
      tune_log_omega = exp(log(tune_log_omega) + pow(i/block,-0.5)*(r_omega-0.234));
      pv_log_omega = pv_log_omega + pow(i/block,-0.5)*(var(log_omega_store.subvec(i-block, i)) - pv_log_omega);
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
    if(i>0 & i%block==0){
      r_alpha = mean(jump_alpha.subvec(i-block, i));
      tune_log_alpha = exp(log(tune_log_alpha) + pow(i/block,-0.5)*(r_alpha-0.234));
      pv_log_alpha = pv_log_alpha + pow(i/block,-0.5)*(var(log_alpha_store.subvec(i-block, i)) - pv_log_alpha);
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
    if(i>0 & i%block==0){
      r_sigma = mean(jump_sigma.subvec(i-block, i));
      tune_log_sigma = exp(log(tune_log_sigma) + pow(i/block,-0.5)*(r_sigma-0.234));
      pv_log_sigma = pv_log_sigma + sqrt(block/i)*(cov(log_sigma_store.rows(i-block, i)) - pv_log_sigma);
      L_log_sigma = chol(pv_log_sigma).t();
    }
    
    // Rcout << "sigma updated" << endl;
    
    // make prediction
    if(i>=burn){
      K_pi_pred = kron(C_pi, H_pred);
      mu_z_pred = X_pred*beta + K_pi_pred*delta_pi;
      sigma2_z_pred = exp(D_pred*log_sigma);
      z_pred = mu_z_pred + sigma2_z_pred%armaNorm(X_pred.n_rows);
      pred_store.row(i-burn) = arma_rpois(exp(z_pred)).t();
    }
    
    prog.increment();
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("z") = z_store.rows(burn, burn+iter-1),
    Rcpp::Named("beta") = beta_store,
    Rcpp::Named("delta_bar") = delta_bar_store, 
    Rcpp::Named("prox") = prox_store,
    Rcpp::Named("kappa_pi") = kappa_pi_store,
    Rcpp::Named("omega")=exp(log_omega_store(span(burn, burn+iter-1))),
    Rcpp::Named("alpha")=exp(log_alpha_store(span(burn, burn+iter-1))),
    Rcpp::Named("sigma")=exp(log_sigma_store.rows(burn,iter+burn-1)),
    Rcpp::Named("pred")=pred_store
  );  
}
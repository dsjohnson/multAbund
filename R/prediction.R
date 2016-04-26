#' @title Predictions for Dirichlet Process Joint Species Distribution Model
#' @description Create MCMC sample from various forms of the predictive distrubution 
#' using the MCMC sample of the parameters
#' @param object A fitted MCMC model object obtained from call to one of the \quote{multAbund}
#' fitting functions.
#' @param pred_list A list similar to data_list in the MCMC fitting functions 
#' except values are used where prediction is desired
#' @param z_pred Alogical indication of whether the posterior of z should 
#' be used (\code{z_pred=FALSE}) or the predictive distribution for z, given the
#' parameter posterior distribution (\code{z_pred=TRUE}; the default). It is ignored
#' for probit occurence models
#' 
#' @export
prediction_sample = function(object, pred_list){
  if(is.null(pred_list$y)){
    num_spec=length(unique(pred_list$data$species))
    Hdelta = apply(fit$delta_bar, 1, 
                   function(x,H,n){kronecker(diag(n),H)%*%x},
                   H = pred_list$H, n=num_spec)
    Xbeta = apply(fit$beta, 1, 
                  function(b,Q){return(Q%*%b)},
                  Q=pred_list$X)
    Sig = apply(X=matrix(log(fit$sigma)), MARGIN=1, 
                function(v,D){return(exp(D%*%v))},
                D=pred_list$D)
    z_new = apply(Sig, 2, function(v){rnorm(length(v), sd=v)})
    lam = exp(Xbeta + Hdelta + z_new)
    pred = matrix(rpois(length(lam), as.vector(lam)), nrow=nrow(lam))
    if(!is.null(object$logit_gamma)){
      gamma = apply(X=object$logit_gamma, MARGIN=1, 
                    function(v,M){return(plogis(M%*%v))},
                    M=pred_list$M)
      zeros = matrix(1-rbinom(length(gamma), size=1, prob=as.vector(gamma)), nrow=nrow(gamma))
      pred = zeros*pred
      return(pred)
    } else return(pred)
  } else{
    prob = pnorm(as.vector(Xbeta + Hdelta))
    pred = matrix(rbinom(length(prob), 1, prob), nrow=nrow(Xbeta))
    return(pred)
  }
}

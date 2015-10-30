## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ------------------------------------------------------------------------
library(multAbund)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(magrittr)
library(reshape2)
library(mvtnorm)

set.seed(555)

## ------------------------------------------------------------------------
### Design
num_spec=20
num_obs=35
num_groups = 5
group_comp = round(num_spec*(c(num_groups:1)/sum(c(num_groups:1))))
env_dat = data.frame(obs=1:num_obs, V1=rnorm(num_obs), V2=rnorm(num_obs), V3=rnorm(num_obs), V4=rnorm(num_obs))
cnt_dat=data.frame(obs=rep(1:num_obs, each=num_spec), species=rep(1:num_spec, num_obs), count=NA)
cnt_dat = merge(cnt_dat, env_dat)
cnt_dat = with(cnt_dat, cnt_dat[order(species, obs), ])

## ------------------------------------------------------------------------
### Model
delta_model = ~V1+V2+V3-1
X_model = ~V1+V2+V3+V4
group = factor(rep(1:num_groups, group_comp))
data_mats = make_data_list(delta_model = delta_model, X_model = X_model, data=cnt_dat, sigma_model=~1)
C_pi = model.matrix(~group-1) 
K_pi = kronecker(C_pi, data_mats$H)
X = data_mats$X

## ------------------------------------------------------------------------
### Parmeters
omega = 1
Omega = omega^2*solve(crossprod(data_mats$H))
delta_pi = rmvnorm(num_groups, sigma=Omega)
delta_pi = as.vector(t(sweep(delta_pi, 2,  apply(delta_pi, 2, mean), "-")))
beta = c(2,1,0,-1,0.5)

z =  X%*%beta + K_pi%*%delta_pi 
cnt_dat$count = rpois(n=length(z), exp(z))

## ------------------------------------------------------------------------
data_mats = make_data_list(counts="count", delta_model = delta_model, 
                           X_model = X_model, data=cnt_dat[cnt_dat$obs<31,], 
                           sigma_model=~1)
pred_mats = make_data_list(counts="count", delta_model = delta_model, 
                           X_model = X_model, data=cnt_dat, 
                           sigma_model=~1)

## ------------------------------------------------------------------------
fit0 = glm(count ~ 1, data=cnt_dat, family="poisson")

#target0=function(n){return(rep(1/n, n))}
alpha_prior = get_opt_alpha_prior(
  num_spec#, target0
)

alpha_ab = alpha_prior$ab
phi_beta = 10
mu_beta = c(fit0$coefficients, rep(0,ncol(X)-1))
phi_omega = 1
df_omega = 1
phi_sigma = 1
df_sigma = 51

# inits=make_inits(
#   data_list = data_mats,
#   phi_beta=10, 
#   mu_beta = c(fit0$coefficients, rep(0,ncol(X)-1)), 
#   phi_omega = 1, 
#   df_omega = 1
# )

inits=sugs(
  data_list = data_mats,
  phi_beta = phi_beta, 
  mu_beta = mu_beta, 
  phi_omega = phi_omega, 
  df_omega = df_omega
)

## ------------------------------------------------------------------------
block=200
burn=20000
iter=100000

fit = mult_abund_pois(
  data_list=data_mats,
  pred_list=pred_mats,
  initial_vals = inits,
  phi_beta = phi_beta, 
  mu_beta = mu_beta, 
  phi_omega = phi_omega, 
  df_omega = df_omega,
  a_alpha=alpha_ab[1],
  b_alpha=alpha_ab[2],
  phi_sigma = phi_sigma,
  df_sigma = df_sigma,
  block = block, 
  begin_group_update=5*block,
  burn = burn, 
  iter = iter
) 

## ---- fig.width=8, fig.height=10-----------------------------------------
# Posterior probability of shared group membership
p1 = ggplot(aes(x=Var1, y=Var2), data=melt(apply(fit$prox, c(1,2), mean))) + geom_tile(aes(fill=value)) + labs(x=NULL, y=NULL) +
  labs(x=NULL, y=NULL) + scale_fill_gradient("", limits=c(0, 1))

# Fitted clustering of species
kappa_dist = table(fit$kappa_pi)
kappa = as.numeric(names(kappa_dist)[kappa_dist==max(kappa_dist)])
d <- as.dist(1-apply(fit$prox, c(1,2), mean))
hc = hclust(d, method="complete") 
dendr    <- dendro_data(hc, type="rectangle") 
clust    <- cutree(hc,k=kappa)                   
clust.df <- data.frame(label=1:20, cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
p2 = ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), size=5) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  theme(
    axis.line.y=element_blank(), axis.text.y=element_blank(), 
    axis.ticks.y=element_blank()
  ) + xlab(NULL)
g1=plot_grid(p1, p2, labels = c("(A) Probability of shared membership", "(B) Estimated clusters"), hjust=0, ncol=1)
print(g1)

# Fit and prediction
ind = cnt_dat$obs<31
pred_data = data.frame(count=cnt_dat$count[!ind], pred=apply(fit$pred, 2, median)[!ind])
p3=ggplot(aes(x=count, y=pred), data=pred_data) + geom_point() + geom_abline(a=0,b=1) + 
  xlab("Unobserved count") + ylab("Predicted count") 
fit_data = data.frame(count=cnt_dat$count[ind], fit=apply(fit$pred, 2, median)[ind])
p4=ggplot(aes(x=count, y=fit), data=fit_data) + geom_point() + geom_abline(a=0,b=1) + 
  xlab("Observed count") + ylab("Fitted count") 
g2=plot_grid(p3, p4, labels = c("(A) Predicted values for unobserved data", "(B) Fitted vlues for observed data"), hjust=0, ncol=1)
print(g2)


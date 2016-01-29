## ----global-options, include=FALSE---------------------------------------
# set global chunk options
library(knitr)
options(digits=2)
options(width=80)
opts_chunk$set(dev="png",
               dev.args=list(type="cairo",antialias="none"),
               dpi=192)

## ------------------------------------------------------------------------
library(multAbund)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(magrittr)
library(reshape2)
library(mvtnorm)
library(viridis)

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
delta_model = ~V1+V2+V3
X_model = ~V1+V2+V3+V4
group = factor(rep(1:num_groups, group_comp))
data_mats = make_data_list(delta_model = delta_model, X_model = X_model, data=cnt_dat, sigma_model=~1)
C_pi = model.matrix(~group-1) 
K_pi = kronecker(C_pi, data_mats$H)
X = data_mats$X

## ------------------------------------------------------------------------
### Parmeters
set.seed(1234)
omega = 1.25
Omega = omega^2*solve(crossprod(data_mats$H))
delta_pi = rmvnorm(num_groups, sigma=Omega)
delta_pi = as.vector(t(sweep(delta_pi, 2,  apply(delta_pi, 2, mean), "-")))
beta = c(2,1,0,-1,0.5)

z =  X%*%beta + K_pi%*%delta_pi 
cnt_dat$count = rpois(n=length(z), exp(z))

## ------------------------------------------------------------------------
data_mats = make_data_list(counts="count", delta_model = delta_model, 
                           X_model = X_model, data=cnt_dat, 
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

prior_parm = list(
  a_alpha=alpha_prior[[1]][1],
  b_alpha=alpha_prior[[1]][2],
  Sigma_beta_inv = crossprod(data_mats$X)/100,
  mu_beta = c(fit0$coefficients, rep(0,ncol(X)-1)),
  phi_omega = 1,
  df_omega = 1,
  phi_sigma = 1,
  df_sigma = 51
)

inits=list(
  beta = with(data_mats, solve(crossprod(X),crossprod(X, log(n+0.5)))),
  groups=c(1:20),
  delta = rep(0, ncol(data_mats$H)*num_spec),
  omega=1,
  sigma = rep(1.0E-4, ncol(data_mats$D)),
  log_alpha=0
)

## ------------------------------------------------------------------------
block=200
burn=10000
iter=50000

fit = mult_abund_pois(
  data_list=data_mats,
  pred_list = data_mats,
  initial_list = inits,
  prior_list = prior_parm,
  block = block, 
  begin_group_update=5*block,
  burn = burn, 
  iter = iter
) 

## ---- fig.width=8, fig.height=10-----------------------------------------
# Posterior probability of shared group membership
p1 = ggplot(aes(x=Var1, y=Var2), data=melt(apply(fit$prox, c(1,2), mean))) + geom_tile(aes(fill=value)) + labs(x=NULL, y=NULL) +
  labs(x=NULL, y=NULL) + scale_fill_viridis("",limits=c(0, 1))

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
fit_data = data.frame(count=cnt_dat$count, fit=apply(fit$pred, 2, mean))
p4=ggplot(aes(x=count, y=fit), data=fit_data) + geom_point() + geom_abline(intercept=0,slope=1) + 
  xlab("Observed count") + ylab("Fitted count") + ggtitle("Fitted model predictions\n")
print(p4)



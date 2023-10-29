library(LaplacesDemon)
## Generate Data:
# Set sigma=0.01, theta0 = 2
# p_i ~ Uniform (0,pi/2), epsilon_{im} ~ N(0, sigma^2), epsilon_{ir} ~ N(0, sigma^2) 
# m_i = sqrt{theta0}*sin(p_i) + epsilon_im, r_i = sqrt{theta0}*sin(p_i) + epsilon_i

sigma= 0.1; theta0=1.5; n= 100
p_true = runif(n,0,pi/2); eps_m = rnorm(n,mean=0,sd=sigma); eps_r = rnorm(n,mean=0,sd=sigma)
m = sqrt(theta0)*sin(p_true) + eps_m; r = sqrt(theta0)*cos(p_true) + eps_r
plot(m,r)

## Model 
# (m_i, r_i)' | p_i, theta, Sigma  ~ N_2( sqrt{theta}*sin(p_i),sqrt{theta}*cos(p_i), Sigma) 
#               p_i ~ Discrete_uniform(seq(0,pi/2, length.out=N=100))
#               theta  ~ Discrete_uniform(seq(theta_L=0,theta_U= 5, length.out=N=100)) # May change this later
#               Sigma ~ Inverse_Wishart(nu = 2, Psi = diag(2))

## Gibbs Sampling
# Global variables
nreps= 500; nu=2; Psi = diag(2); N=200; theta_L =0; theta_U=5
grid_theta = seq(theta_L, theta_U, length.out=N)
grid_p=seq(0, pi/2, length.out=N)
# Unknown parameters: theta, p=(p_1,..., p_n), Sigma
# Storage variables:
theta_store = numeric(nreps); p_store = matrix(0,nrow=nreps,ncol=n); 
Sigma_store = array(0,dim = c(2,2,nreps))
# Initialization: update theta first, so initialize p  and Sigma (draws from prior)
p=sample(seq(0,pi/2, length.out=N), n, replace=T); Sigma = rinvwishart(nu, Psi)

for (iter in 1:nreps)
{
  iSigma = solve(Sigma)
  # Update theta
  theta_prob_log= mapply(function(theta){res = rbind(m - sqrt(theta)*sin(p),r - sqrt(theta)*cos(p)); 
  -0.5*sum(diag(t(res)%*%iSigma%*%res))},grid_theta)
  unnorm_prob = exp( theta_prob_log - max(theta_prob_log))
  norm_prob = unnorm_prob / sum(unnorm_prob)
  theta =  sample(grid_theta, 1, prob = norm_prob)
  theta_store[iter] = theta
  
  # Update p
  for (i in 1:n)
  {
    p_prob_log= mapply(function(x){res = c(m[i] - sqrt(theta)*sin(x),r[i] - sqrt(theta)*cos(x)); 
    -0.5*res%*%iSigma%*%res},grid_p)
    unnorm_prob = exp( p_prob_log - max(p_prob_log))
    norm_prob = unnorm_prob / sum(unnorm_prob)
    p_store[iter,i] = sample(grid_p, 1, prob = norm_prob)
  }
  p = p_store[iter,]
  
  # Update Sigma
  res = rbind(m - sqrt(theta)*sin(p),r - sqrt(theta)*cos(p));  S <- res%*%t(res)
  # Updated parameters for inverse Wishart
  nu_star =nu + n; Psi_star <- Psi + S 
  Sigma = rinvwishart(nu_star, Psi_star)
  Sigma_store[,,iter] = Sigma
  
}

hist(theta_store,50)
plot(apply(p_store,2,mean),p_true)
abline(coef = c(0,1))

length(m) ; length(r)
N = 200 ; n = 100
grid_p <- seq(0,pi/2, length.out=N) ; p=sample(grid_p, 100, replace=T)
theta = 2
mean <- sqrt(theta)*sin(p) ; mean ; length(mean)
res = rbind(m - mean,r - mean) ; res ; dim(res) ;  iSigma = solve(Sigma)
sum(diag(t(res)%*%iSigma%*%res))

library(LaplacesDemon)
## Generate Data:
   # Set sigma=0.01, theta0 = 2
   # p_i ~ Uniform (0,pi/2), epsilon_{im} ~ N(0, sigma^2), epsilon_{ir} ~ N(0, sigma^2) 
   # m_i = sqrt{theta0}*sin(p_i) + epsilon_im, r_i = sqrt{theta0}*sin(p_i) + epsilon_i

    sigma= 0.1; theta0=1.5; n= 100
    p_true = runif(n,0,pi/2); eps_m = rnorm(n,mean=0,sd=sigma); eps_r = rnorm(n,mean=0,sd=sigma)
    m = sqrt(theta0)*sin(p_true) + eps_m; r = sqrt(theta0)*cos(p_true) + eps_r
    plot(m,r)
    
## Model 
    # (m_i, r_i)' | p_i, theta, Sigma  ~ N_2( sqrt{theta}*sin(p_i),sqrt{theta}*cos(p_i), Sigma) 
    #               p_i ~ Discrete_uniform(seq(0,pi/2, length.out=N=100))
    #               theta  ~ Discrete_uniform(seq(theta_L=0,theta_U= 5, length.out=N=100)) # May change this later
    #               Sigma ~ Inverse_Wishart(nu = 2, Psi = diag(2))

## Gibbs Sampling
    # Global variables
    nreps= 500; nu=2; Psi = diag(2); N=200; theta_L =0; theta_U=5
    grid_theta = seq(theta_L, theta_U, length.out=N)
    grid_p=seq(0, pi/2, length.out=N)
    # Unknown parameters: theta, p=(p_1,..., p_n), Sigma
    # Storage variables:
    theta_store = numeric(nreps); p_store = matrix(0,nrow=nreps,ncol=n); 
    Sigma_store = array(0,dim = c(2,2,nreps))
    # Initialization: update theta first, so initialize p  and Sigma (draws from prior)
    p=sample(seq(0,pi/2, length.out=N), n, replace=T); Sigma = rinvwishart(nu, Psi)
    
    for (iter in 1:nreps)
    {
      iSigma = solve(Sigma)
      # Update theta
        theta_prob_log= mapply(function(theta){res = rbind(m - sqrt(theta)*sin(p),r - sqrt(theta)*cos(p)); 
         -0.5*sum(diag(t(res)%*%iSigma%*%res))},grid_theta)
       unnorm_prob = exp( theta_prob_log - max(theta_prob_log))
       norm_prob = unnorm_prob / sum(unnorm_prob)
       theta =  sample(grid_theta, 1, prob = norm_prob)
       theta_store[iter] = theta
      
       # Update p
       for (i in 1:n)
       {
         p_prob_log= mapply(function(x){res = c(m[i] - sqrt(theta)*sin(x),r[i] - sqrt(theta)*cos(x)); 
         -0.5*res%*%iSigma%*%res},grid_p)
         unnorm_prob = exp( p_prob_log - max(p_prob_log))
         norm_prob = unnorm_prob / sum(unnorm_prob)
         p_store[iter,i] = sample(grid_p, 1, prob = norm_prob)
       }
       p = p_store[iter,]
       
       # Update Sigma
        res = rbind(m - sqrt(theta)*sin(p),r - sqrt(theta)*cos(p));  S <- res%*%t(res)
        # Updated parameters for inverse Wishart
        nu_star =nu + n; Psi_star <- Psi + S 
        Sigma = rinvwishart(nu_star, Psi_star)
        Sigma_store[,,iter] = Sigma
      
    }
  
  hist(theta_store,50)
  plot(apply(p_store,2,mean),p_true)
  abline(coef = c(0,1))
 
  length(m) ; length(r)
  N = 200 ; n = 100
  grid_p <- seq(0,pi/2, length.out=N) ; p=sample(grid_p, 100, replace=T)
  theta = 2
  mean <- sqrt(theta)*sin(p) ; mean ; length(mean)
  res = rbind(m - mean,r - mean) ; res ; dim(res) ;  iSigma = solve(Sigma)
  sum(diag(t(res)%*%iSigma%*%res))
  

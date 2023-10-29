#Two dimensional theta simulation:

library(LaplacesDemon)

##Generate Data:
sigma = 0.1; theta1_0 = 4; theta2_0 = 1; n = 100; p_true = runif(n, 0,pi/2)
m = sqrt(theta1_0)*cos(p_true) +  rnorm(n, mean=0, sd=sigma)
r = sqrt(theta2_0)*sin(p_true) + rnorm(n, mean=0, sd=sigma)
plot(m, r)

##Gibbs-Sampling
#Global variables
nreps = 10000; nu = 2; Psi = diag(2); N = 200
grid_theta = seq(0, 5, length.out=N)
grid_p = seq(0, 2*pi, length.out=N)
#storage variables          
theta1_store = numeric(nreps); theta2_store = numeric(nreps)
p_store = matrix(0, nrow=nreps, ncol=n) ; Sigma_store = array(0, dim=c(2, 2, nreps))
#Initialization
p = sample(grid_p, n, replace=T); Sigma = rinvwishart(nu, Psi); theta2 = runif(1, min = 0, max = 5)

for (iter in 1:nreps) {
  iSigma <- solve(Sigma); L <- chol(Sigma) ; ihalfSigma =  solve(t(L)) %*% solve(L)
  # Update theta1
  theta1_prob_log = mapply(function(theta1) {
    res = rbind(m - sqrt(theta1)*cos(p), r - sqrt(theta2)*sin(p))
    -0.5 * sum(rowSums((t(res) %*% ihalfSigma)^2))
  }, grid_theta)
  unnorm_prob1 = exp(theta1_prob_log - max(theta1_prob_log))
  norm_prob1 = unnorm_prob1 / sum(unnorm_prob1)
  theta1 = sample(grid_theta, 1, prob=norm_prob1)
  theta1_store[iter] = theta1
  
  # Update theta2
  theta2_prob_log = mapply(function(theta2) {
    res = rbind(m - sqrt(theta1)*cos(p), r - sqrt(theta2)*sin(p))
    -0.5 * sum(rowSums((t(res) %*% ihalfSigma)^2))
  }, grid_theta)
  unnorm_prob2 = exp(theta2_prob_log - max(theta2_prob_log))
  norm_prob2 = unnorm_prob2 / sum(unnorm_prob2)
  theta2 = sample(grid_theta, 1, prob=norm_prob2)
  theta2_store[iter] = theta2
  
  # Update p
  for (i in 1:n) {
    p_prob_log = mapply(function(x) {
      res = c(m[i] - sqrt(theta1)*cos(x), r[i] - sqrt(theta2)*sin(x))
      -0.5 * res %*% iSigma %*% res
    }, grid_p)
    unnorm_prob = exp(p_prob_log - max(p_prob_log))
    norm_prob = unnorm_prob / sum(unnorm_prob)
    p_store[iter, i] = sample(grid_p, 1, prob=norm_prob)
  }
  p = p_store[iter, ]
  
  # Update Sigma
  res = rbind(m - sqrt(theta1)*cos(p), r - sqrt(theta2)*sin(p))
  S = res %*% t(res)
  nu_star = nu + n
  Psi_star = Psi + S
  Sigma = rinvwishart(nu_star, Psi_star)
  Sigma_store[, , iter] = Sigma
}

##Results Visualization
hist(theta1_store, 50, main="Histogram of Theta1")
hist(theta2_store, 50, main="Histogram of Theta2")
plot(apply(p_store, 2, mean), p_true)
abline(coef=c(0, 1))

# Trace plot for theta1
plot(theta1_store, type="l", col="blue", main="Trace Plot for Theta1", ylab="Theta1", xlab="Iteration")
abline(h = true_theta1, col="red", lwd=2)

# Trace plot for theta2
plot(theta2_store, type="l", col="darkgreen", main="Trace Plot for Theta2", ylab="Theta2", xlab="Iteration")
abline(h = true_theta2, col="red", lwd=2)


#----Microbenchmark:
library(microbenchmark)

sum_diagonal_approach = function() { mapply(function(theta1) {
  res = rbind(m - sqrt(theta1)*cos(p), r - sqrt(theta2)*sin(p))
  -0.5 * sum(diag(t(res) %*% solve(Sigma) %*% res))
}, grid_theta)}

euclidean_approach = function() { mapply(function(theta1) {
  res = rbind(m - sqrt(theta1)*cos(p), r - sqrt(theta2)*sin(p))
  -0.5 * sum(rowSums((t(res) %*% (solve(t(chol(Sigma))) %*% solve(chol(Sigma))))^2))
}, grid_theta)}

microbenchmark(sum_diagonal_approach, euclidean_approach, times = 1000)


m = sqrt(theta1_0)*cos(p_true) 
r = sqrt(theta2_0)*sin(p_true) 
plot(m, r)

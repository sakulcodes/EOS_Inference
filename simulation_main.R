# Libraries
library(MASS) # for mvrnorm function
library(mvtnorm)
library(tmvtnorm)
library(LaplacesDemon)
library(ggplot2)
library(MCMCpack)
library(base)

# ---------------------------------------------DATA GENERATION ----------------------------------------------------------
# Given theta, generate_data function generates (M, R) data
generate_data <- function(theta, n_points = 100) {
  pi_val <- runif(n_points, 0, pi / 2) # Uniform Sampling #narrow it down (0, pi/4)
  M_true <- sqrt(theta) * cos(pi_val)
  R_true <- sqrt(theta) * sin(pi_val) 
  return(data.frame(M_true, R_true, pi_val))
}

# Simulate Data with Gaussian scatter
simulate_observations <- function(true_data, Sigma) {
  n <- nrow(true_data)
  observed_data <- data.frame(M = rep(0, n), R = rep(0, n))
  for(i in 1:n) {
    error <- mvrnorm(n = 1, mu = c(0, 0), Sigma = Sigma)
    observed_data$M[i] <- true_data$M_true[i] + error[1]
    observed_data$R[i] <- true_data$R_true[i] + error[2]
  }
  return(observed_data)
}

# True Values
theta_true <- 1 # True value of theta
Sigma_true <- matrix(c(0.001, 0, 0, 0.001), 2, 2) 

# Simulate data based on known theta_true and Sigma_true
true_data <- generate_data(theta_true)
observed_data <- simulate_observations(true_data, Sigma_true)
plot(observed_data)
true_p_values <- true_data$pi_val


#pi_val, should be defined in a grid from 0 to 2pi, cannot use observations from data. 
#M_expected and R_expected functions
M_expected <- function(pressure_values, theta) {
  return(sqrt(theta) * cos(pressure_values))
}

R_expected <- function(pressure_values, theta) {
  return(sqrt(theta) * sin(pressure_values))
}

#Initialize theta and p_values
#theta <- theta_true
sigma <- Sigma_true

#Define the discerete grid for pi values
L <- 0          
U <- pi/2     
Na <- 99 #no. of intervals
delta_p <- (U - L) / Na #step size
pi_values <- L + seq(0, Na) * delta_p #grid for theta values

#Define the discrete grid for theta values
Lower_limit <- 0          
Upper_limit <- 5
No <- 99 #no. of intervals
delta_theta <- (Upper_limit - Lower_limit) / No #step size
theta_values <- Lower_limit + seq(0, No) * delta_theta #grid for theta values


log_likelihood <- function(data, theta, Sigma) {
  n <- nrow(data)
  log_likelihood_value <- 0  
  for(i in 1:n) {
    pi_val <- pi_values[i]
    M_val <- M_expected(pi_val, theta)
    R_val <- R_expected(pi_val, theta)
    mu <- c(M_val, R_val)
    obs <- c(data$M[i], data$R[i])
    log_likelihood_value <- log_likelihood_value + dmvnorm(x = obs, mean = mu, sigma = Sigma, log = TRUE)
  }
  return(log_likelihood_value)
}

print(log_likelihood(observed_data, theta_values[1], Sigma_true))


#Initialize parameters for sigma 
nu = 2 # minimal non-informative choice
Psi = diag(2) #2*2 identity matrix

#--------------------------------------Sampling Steps --------------------------------------------

#---------For theta

log_prior_theta <- function(theta, Lower_limit, Upper_limit) {
    return(-log(Upper_limit - Lower_limit))
}
print(log_prior_theta(theta_values, Lower_limit, Upper_limit))
      
      
sample_theta_directly <- function(data, Sigma, theta_values, Lower_limit, Upper_limit) {
  # Compute unnormalized posterior for each theta value
  unnormalized_log_posteriors <- sapply(theta_values, function(theta_val) {
    log_likelihood_val = log_likelihood(data, theta_val, Sigma)
    log_prior_val = log_prior_theta(theta_val, Lower_limit, Upper_limit)
    return(log_likelihood_val + log_prior_val)
  })
  
  # Convert log values to probabilities using the log-sum-exp trick
  max_log = max(unnormalized_log_posteriors)
  unnormalized_posteriors = exp(unnormalized_log_posteriors - max_log)
  normalization_constant = sum(unnormalized_posteriors)
  normalized_posteriors = unnormalized_posteriors / normalization_constant
  
  #Sample from discrete distribution
  sampled_index = sample(1:length(theta_values), 1, prob = normalized_posteriors)
  
  return(theta_values[sampled_index])
}


#-----------For Sigma
sample_sigma <- function(data, pressure_values, theta, nu, Psi) {
  N = nrow(data)
  residuals <- matrix(0, N, 2)
  for (i in 1:N) {
    residuals[i, ] <- c(data$M[i], data$R[i]) - c(M_expected(pressure_values[i], theta), R_expected(pressure_values[i], theta))
  }
  S <- t(residuals) %*% residuals 
  
  # Updated parameters for inverse Wishart
  nu_star <- nu + N
  Psi_star <- Psi + S
  
  # Sample Sigma from the inverse Wishart distribution
  return(rinvwishart(nu_star, Psi_star))
}

print(sample_sigma(data = observed_data, pressure_values = pi_values, theta = theta_true, nu, Psi))

#----------For P_i

log_prior_pi <- function(pressure_values, L, U) {
  return(-log(U-L))
}

sample_pi_directly <- function(data, theta, sigma, L, U, pressure_values) {
  # Compute unnormalized log posteriors for each pi value
  unnormalized_posteriors <- sapply(pressure_values, function(pressure_values) {
    log_likelihood_val = log_likelihood(data, theta, sigma) 
    log_prior_val = log_prior_pi(pressure_values, L, U)
    return(log_likelihood_val + log_prior_val)
  })
  
  # Convert unnormalized log posteriors to probabilities
  unnormalized_probs = exp(unnormalized_posteriors - max(unnormalized_posteriors))
  normalized_probs = unnormalized_probs / sum(unnormalized_probs)
  
  # Sample from the discrete distribution of pi_values
  sampled_pi = sample(pressure_values, 1, prob = normalized_probs)
  
  return(sampled_pi)
}

print(sample_pi_directly(observed_data, theta, sigma, L, U, pi_values))


#------------------------------------GIBBS SAMPLING STEPS --------------------------------------------
# Number of Gibbs iterations 
n_iterations <- 1000
N = nrow(observed_data) #total number of data points

# Initialize variables to store results
theta_estimates <- numeric(n_iterations)
Sigma_estimates <- list(n_iterations)
p_i_estimates <- list(n_iterations)

#Intialize the parameters value:
theta_estimates[1] = 1; Sigma_estimates[[1]] = Sigma_true; p_i_estimates[[1]] = pi_values

for (iteration in 2:n_iterations) {
  
  theta_estimates[iteration] = sample_theta_directly(theta = theta_values, data = observed_data, Sigma = Sigma_estimates[[iteration-1]], Lower_limit = Lower_limit, Upper_limit = Upper_limit)
  Sigma_estimates[[iteration]] = sample_sigma(data = observed_data, pressure_values = p_i_estimates[[iteration-1]], theta = theta_estimates[iteration], nu, Psi)
  p_i_estimates[[iteration]] = numeric(N)
  for (i in 1:N) {
    p_i_estimates[[iteration]][i] = sample_pi_directly(data = observed_data, theta = theta_estimates[iteration], sigma = Sigma_estimates[[iteration]], L, U, pressure_values =  pi_values)
  }
}

# Burn-in and thinning
burn_in = 100  # for example
thin = 10  # for example

theta_samples <- theta_estimates[-seq(1, burn_in)]
Sigma_samples <- Sigma_estimates[-seq(1, burn_in)]
p_i_samples < p_i_estimates[-seq(1, burn_in)]

theta_samples <- theta_samples[seq(1, length(theta_samples), by = thin)]
Sigma_samples <- Sigma_samples[seq(1, length(Sigma_samples), by = thin)]
p_i_samples <- p_i_samples[seq(1, length(p_i_samples), by = thin)]




#------------------------------------------ANALYZING RESULT------------------------------------------------------------------------
#For theta:
#mean and standard deviation 
mean_theta_estimates <- mean(theta_estimates)
sd_theta_estimates <- sd(theta_estimates)
cat("Mean of estimated theta:", mean_theta_estimates, "\n")
cat("Standard Deviation of estimated theta:", sd_theta_estimates, "\n")


#density plot
theta_df <- data.frame(theta_estimates)
ggplot(theta_df, aes(x = theta_estimates)) +
  geom_density(fill = "blue", alpha = 0.6) +
  geom_vline(xintercept = sqrt(theta_true), color = "red", linetype = "dashed") +
  labs(title = "Density Plot of Estimated Theta", x = "Theta") +
  theme_minimal() + scale_x_continuous(limits = c(-4, 4))

 #95% Confidence Interval for theta
ci <- quantile(theta_estimates, c(0.025, 0.975)) 
print(ci)
#For Sigma:
print(Sigma_estimates)
print(cov(observed_data))


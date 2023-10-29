library(LaplacesDemon)
library(MASS)

# Generate Data:
# obs_data <- (mi_obs, ri_obs) ~  {m(p_i_true, theta_true),r(p_i_true, theta_true)} +  true_error
# true_error ~ N2(0,Sigma_true)
# Set Sigma_true = cov(obs_data)

# ----------------------observed - Riley and Miller (NICER DATASET)
#(M,R) data-set
riley_data <- read.delim("/Users/sakul/Desktop/RESEARCH/PRL_rev/riley/ST_PST/run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1post_equal_weights.dat", sep=" ", header=FALSE)[, c(13, 9)]
miller_data <- read.table("/Users/sakul/Desktop/RESEARCH/PRL_rev/miller/J0030_3spot_RM.txt", header=FALSE)[sample(1:10000, 12242, replace=TRUE),]
combined_data <- rbind(na.omit(setNames(riley_data, c("R", "M"))), setNames(miller_data, c("R", "M")))
observed_M <- combined_data$M ; observed_R <- combined_data$R
obs_cov_mat <- cov(combined_data) #covariance matrix

#generate data:
n = length(observed_M)
correlated_noise = mvrnorm(n, mu = c(0, 0), Sigma = obs_cov_mat)
m_observed = observed_M + correlated_noise[,1] ; r_observed = observed_R + correlated_noise[,2]
plot(m_observed, r_observed)

#model 
# (m_i, r_i)' | p_i, theta, Sigma  ~ N_2( {m(p_i, theta),r(p_i, theta)} , Sigma) 
#               p_i ~ uniform(pi-prior-data)
#               theta  ~ Uniform(theta-prior-data) 
#               Sigma ~ Inverse_Wishart(nu = 2, Psi = diag(2))

#grid for theta 
prior <- read.table("/Users/sakul/Desktop/ldmcrust/data/eft_pnm32_kfr16_par_all.dat")[, -c(5, 6)]
grid_theta <- prior ; grid_theta <- grid_theta[1:100,]  #for psuedo purposes, sampling only first 100 values as orcale also has only data points.
log_prior_theta = -log(nrow(grid_theta))

# ------------------------ {M(pi,theta), R(pi,theta)} - implicit relationship - RESULT MATRIX
#The FORTRAN oracle should solve (M(pi,theta), R(pi, theta)) FOR EVERY COMBINATION in the grid values for theta and p

#Creating an oracle matrix of dimensions = k (number of theta) * n(number of central densities / central pressure) for every combination in the grid. 
max_rows <- 0 ; num_files <- 100 ; base_path <- "/Users/sakul/Desktop/ldmcrust/res_eos_pnm32/eft_pnm32_"
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)
  max_rows <- max(max_rows, nrow(data_eos))
}
M_values <- matrix(NA, nrow = max_rows, ncol = num_files) ; R_values <- matrix(NA, nrow = max_rows, ncol = num_files)

# Fill the matrices
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)
  M_values[1:nrow(data_eos), i] <- data_eos$V2
  R_values[1:nrow(data_eos), i] <- data_eos$V3
}

M_values; R_values #oracle matrix formulation, each row corresponds to a central density and each column is the solution of a theta value. 
dim(M_values) ; M_values[300, 2] #should give NA
#Note - FORTRAN EOS solution needs to be extended if we decide to extend the grid. 

#Replacing the NAs values in every column with the mean of Non-NA values in that column (mean imputation) - just for now, will think about this later.
M_values <- apply(M_values, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  return(col)
})

R_values <- apply(R_values, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  return(col)
})

#grid for central density / central pressure
all_densities <- c()
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)
  all_densities <- c(all_densities, data_eos$V1)
}
unique_densities <- sort(unique(all_densities)) ; unique_densities
grid_p <- unique_densities ; log_prior_p = -log(length(grid_p))


# -------------------------------Gibbs-Sampling-Setup:
nreps = 500; nu = 2; Psi = diag(2) ; n = length(observed_M)
#storage variables:
theta_store = matrix(0, nrow = nreps, ncol = 8) ; p_store = matrix(0, nrow = nreps, ncol = n) ; Sigma_store = array(0, dim = c(2, 2, nreps))
p = sample(grid_p, n, replace = TRUE) ; Sigma = rinvwishart(nu, Psi)  #initialization   

for (iter in 1:nreps) {
  iSigma = solve(Sigma)
  # Update theta vector
  theta_prob_log = sapply(1:nrow(grid_theta), function(i) {  # Loop over columns (theta values)
    expected_m = M_values[, i] ; expected_r = R_values[, i]
    res = rbind(m_observed - expected_m[p], r_observed - expected_r[p])
    -0.5 * sum(diag(t(res) %*% iSigma %*% res))
  })
  unnorm_prob = exp(theta_prob_log - max(theta_prob_log))
  norm_prob = unnorm_prob / sum(unnorm_prob)
  sampled_row_idx = sample(1:nrow(grid_theta), 1, prob = norm_prob)
  sampled_theta = grid_theta[sampled_row_idx, ]
  theta_store[iter, ] = sampled_theta
  
  # Update p
  theta_idx = which(apply(grid_theta, 1, function(row) all(row == theta)))
  
  for (i in 1:n)
  {
    p_prob_log = mapply(function(p_val) {
      # Get the expected m and r values for the current p and previously sampled theta
      expected_m = M_values[p_val, theta_idx]
      expected_r = R_values[p_val, theta_idx]
      res = c(m_observed[i] - expected_m, r_observed[i] - expected_r)
      -0.5 * res %*% iSigma %*% res + log_prior_p
    }, grid_p)
    unnorm_prob = exp(p_prob_log - max(p_prob_log))
    norm_prob = unnorm_prob / sum(unnorm_prob)
    p_store[iter, i] = sample(grid_p, 1, prob = norm_prob)
  }
  p = p_store[iter,]
  
  #update sigma 
  #find residuals based on previously sampled theta and p
  residuals_m = m_observed - M_values[p, theta_idx] ; residuals_r = r_observed - R_values[p, theta_idx] 
  res = rbind(residuals_m, residuals_r)
  S = res %*% t(res)
  # Updated the parameters for inverse wishart
  nu_star = nu + n
  Psi_star = Psi + S
  Sigma = rinvwishart(nu_star, Psi_star)
  Sigma_store[, , iter] = Sigma
}

#Trace plots for theta:
plot(theta_store[, 1], type="l", col="blue", xlab="Iteration", ylab="Theta[1]") #Trace for first parameter of theta:
par(mfrow=c(4,2))  #Trace for all 8 parameters
for (i in 1:8) {
  plot(theta_store[, i], type="l", col="blue", xlab="Iteration", ylab=paste("Theta[", i, "]", sep=""))
}

#Do some more goodness of fit checks:



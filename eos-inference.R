library(LaplacesDemon)
library(MASS) 
library(ggplot2)

# Generate Data:
# obs_data <- (mi_obs, ri_obs) ~  {m(p_i_true, theta_true),r(p_i_true, theta_true)} +  true_error
# true_error ~ N2(0,Sigma_true)
# Set Sigma_true = cov(obs_data)

# ---------------------observed - Riley and Miller (NICER DATASET)
#(M,R) data-set
riley_data <- read.delim("/Users/sakul/Desktop/RESEARCH/PRL_rev/riley/ST_PST/run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1post_equal_weights.dat", sep=" ", header=FALSE)[, c(13, 9)]
miller_data <- read.table("/Users/sakul/Desktop/RESEARCH/PRL_rev/miller/J0030_3spot_RM.txt", header=FALSE)[sample(1:10000, 12242, replace=TRUE),]
combined_data <- rbind(na.omit(setNames(riley_data, c("R", "M"))), setNames(miller_data, c("R", "M")))
observed_M <- combined_data$M ; observed_R <- combined_data$R
correlation <- cor(observed_M, observed_R); print(correlation)
plot(observed_R, observed_M, xlab = "Observed Radius", ylab = "Observed Mass", main ="Observed Data")
obs_cov_mat <- cov(combined_data) ; obs_cov_mat#covariance matrix

#generate data:
m_observed = observed_M ; r_observed = observed_R 
plot(r_observed, m_observed, main ="Observed Data")

#Bayesian Parametric Model 
# (m_i, r_i)' | p_i, theta, Sigma  ~ N_2( {m(p_i, theta),r(p_i, theta)} , Sigma) 
#              theta  ~ Discrete Uniform(theta_1,....,theta_K) 
#               p_i | theta ~ Discrete Uniform {Pmin(theta),....,Pmax(theta)} independently for i = 1 to n
#               Sigma ~ Inverse_Wishart(nu = 2, Psi = diag(2))


# ------------------------ ORACLE MATRIX  - {M(pi,theta), R(pi,theta)} - Implicit Relationship - Solves the EOS & TOV:
#The FORTRAN oracle should solve (M(pi,theta), R(pi, theta)) FOR EVERY COMBINATION in the grid values for theta and p. 
#Imagine a N by K matrix, N = no. of central densities/central pressure, K = no.of theta 

#creating a lilst of all_data to work with -> think of this as the oracle matrix
num_files <- 500 # Assuming you have 100 files
base_path <- "/Users/sakul/Desktop/ldmcrust/res_eos_pnm32/eft_pnm32_"
all_data <- vector("list", num_files)
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)[, 1:3]  # Read only the first three columns
  all_data[[i]] <- data_eos
}

# -----> Parametric Solution Curve visualization for 1st, 2nd, 3rd and 100th EOS overlayed

plot( all_data[[1]][,3], all_data[[1]][,2], xlab = "Radius", ylab = "Mass", col = "blue") ; points(all_data[[2]][,3], all_data[[2]][,2], col = "red")
points(all_data[[3]][,3], all_data[[3]][,2], col = "green")
points(all_data[[100]][,3], all_data[[100]][,2], col = "pink")


#check what linear combinations result in this ({ai},{bi}) grid. - WORK ON THIS:
#grid for theta: #
prior <- read.table("data/eft_pnm32_kfr16_par_all.dat")[, -c(5, 6)]
grid_theta <- prior ; grid_theta <- grid_theta[1:500,]  
#for pseudo purposees - taking grid for only the first 100 values as only 100 data points are loaded in the oracle matrix. 


#the firs column of nsmr.dat file is the central densities
#not using this code:
# #general grid for central density(nc)
# all_densities <- c()
# for (i in 1:num_files) {
#   file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
#   data_eos <- read.table(file_name)
#   all_densities <- c(all_densities, data_eos$V1)
# }
# unique_densities <- sort(unique(all_densities), decreasing = TRUE) ; unique_densities
# grid_p <- unique_densities

# ------------------------------Gibbs-Sampling-Setup:
nreps = 200; nu = 2; Psi = diag(2) ; n = length(m_observed)
#storage variables:
theta_store = matrix(0, nrow = nreps, ncol = 8) ; p_store = matrix(0, nrow = nreps, ncol = n) ; Sigma_store = array(0, dim = c(2, 2, nreps))
#Initialize theta (from grid_theta take the row closest to the mean)
theta_grid <- as.matrix(grid_theta)
mean = colMeans(theta_grid) ; distances = apply(theta_grid, 1, function(row) sum((row - mean)^2))  # Calculate squared distances
theta = theta_grid[which.min(distances), ]

#Initialize Sigma:
Sigma = rinvwishart(nu, Psi)  #initialization

for (iter in 1:nreps) {
  print(iter)
  iSigma = solve(Sigma)
  iter = 2
  #Update p
  theta_idx = which(apply(theta_grid, 1, function(row) all(row == theta)))
  for (i in 1:n) {
    p_prob_log = mapply(function(p_val) {
      row_index = which(all_data[[theta_idx]][,1] == p_val)
      expected_m = all_data[[theta_idx]][row_index,2]
      expected_r = all_data[[theta_idx]][row_index,3]
      res = c(m_observed[i] - expected_m, r_observed[i] - expected_r)
      return(-0.5 * res %*% iSigma %*% res)
    }, all_data[[theta_idx]][,1])
    unnorm_prob = exp(p_prob_log - max(p_prob_log))
    norm_prob = unnorm_prob / sum(unnorm_prob)
    p_store[iter, i] = sample(all_data[[theta_idx]][,1], 1, prob = norm_prob)
  }
  p = p_store[iter,]
  
  #update theta:
  theta_prob_log = mapply(function(idx) {
    row_indices = match(unique(p), all_data[[idx]][,1])
    valid_idx = !is.na(row_indices)
    expected_m = all_data[[idx]][row_indices[valid_idx],2]
    expected_r = all_data[[idx]][row_indices[valid_idx],3]
    print(idx)
    res = rbind(m_observed - expected_m, r_observed - expected_r)
    -0.5 * sum(diag(t(res) %*% iSigma %*% res))
  }, seq_len(length(all_data)))
  
  unnorm_prob = exp(theta_prob_log - max(theta_prob_log))
  norm_prob = unnorm_prob / sum(unnorm_prob)
  sampled_theta = theta_grid[sample(1:nrow(grid_theta), 1, prob = norm_prob), ]
  theta_store[iter, ] = sampled_theta
  theta = sampled_theta
  
  #update sigma:
  #find residuals based on previously sampled theta and p
  theta_idx = which(apply(theta_grid, 1, function(row) all(row == sampled_theta)))
  row_indexx = match(unique(p),all_data[[theta_idx]][,1])
  valid_idx = !is.na(row_indexx)
  residuals_m = m_observed - all_data[[theta_idx]][row_indexx[valid_idx], 2]
  residuals_r = r_observed - all_data[[theta_idx]][row_indexx[valid_idx], 3]
  res = rbind(residuals_m, residuals_r)
  S = res %*% t(res)
  # Updated the parameters for inverse wishart
  nu_star = nu + n
  Psi_star = Psi + S
  Sigma = rinvwishart(nu_star, Psi_star)
  Sigma_store[, , iter] = Sigma
}

# Assuming theta_store is a matrix with rows as iterations and columns as theta parameters
num_iters <- 2

#visualizing trace for ai coefficients/parameters:
plot(theta_store[, 1], type = "l", col = 1,
     xlab = "Iteration", ylab = "Theta Value", 
     main = "Trace Plot of a0")

plot(theta_store[, 2], type = "l", col = 1,
     xlab = "Iteration", ylab = "Theta Value", 
     main = "Trace Plot of a1")

plot(theta_store[, 3], type = "l", col = 1,
     xlab = "Iteration", ylab = "Theta Value", 
     main = "Trace Plot of a2")

plot(theta_store[, 4], type = "l", col = 1,
     xlab = "Iteration", ylab = "Theta Value", 
     main = "Trace Plot of a3")



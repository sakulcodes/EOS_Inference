library(LaplacesDemon)
library(MASS) 
library(ggplot2)

#Generate Data:
# ---------------------observed - Riley and Miller (NICER DATASET)
#(M,R) data-set
riley_data <- read.delim("data/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1post_equal_weights.dat", sep=" ", header=FALSE)[, c(13, 9)]
miller_data <- read.table("data/J0030_3spot_RM.txt", header=FALSE)[sample(1:10000, 12242, replace=TRUE),]
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

#creating a list of all_data -> oracle matrix:
num_files <- 500 #no.of files to process
base_path <- "data/res_eos_pnm32/eft_pnm32_"
all_data <- vector("list", num_files)

#loading the 500 to 1000 theta entries:
for (i in 500:(500 + num_files - 1)) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)[, 1:3]  # Read only the first three columns
  all_data[[i - 499]] <- data_eos  # Adjust index for all_data list
}

# -----> Example Visualizaiton: Parametric Solution Curve visualization for 1st, 2nd, 3rd and 100th EOS overlayed (not needed for computation)

plot( all_data[[1]][,3], all_data[[1]][,2], xlab = "Radius", ylab = "Mass", col = "blue") ; points(all_data[[2]][,3], all_data[[2]][,2], col = "red")
points(all_data[[3]][,3], all_data[[3]][,2], col = "green")
points(all_data[[100]][,3], all_data[[100]][,2], col = "pink")


#grid for theta: 
prior <- read.table("data/eft_pnm32_kfr16_par_all.dat")[, -c(5, 6)] #the 5th and 6th columns are zero. 
grid_theta <- as.matrix(prior)
grid_theta <- grid_theta[500:1000,]  


# ------------------------------Gibbs-Sampling-Setup:
nreps = 200; nu = 2; Psi = diag(2) ; n = length(m_observed)
#storage variables:
theta_store = matrix(0, nrow = nreps, ncol = 8) ; p_store = matrix(0, nrow = nreps, ncol = n) ; Sigma_store = array(0, dim = c(2, 2, nreps))

#Initialize theta (from grid_theta take the row closet to the mean)
theta_grid <- as.matrix(grid_theta)
mean = colMeans(theta_grid) ; distances = apply(theta_grid, 1, function(row) sum((row - mean)^2))  # Calculate euclidean distance for each theta from the mean.
theta = theta_grid[which.min(distances), ] #initialized theta

#Initialize Sigma:
Sigma = rinvwishart(nu, Psi)  #initialization

for (iter in 1:nreps) {
  
  print(iter)
  iSigma = solve(Sigma) ; L <- chol(Sigma) ; ihalfSigma =  solve(t(L)) %*% solve(L)
  
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
    row_indices = match(p, all_data[[idx]][,1])
    valid_idx = !is.na(row_indices)
    expected_m = all_data[[idx]][row_indices[valid_idx],2]
    expected_r = all_data[[idx]][row_indices[valid_idx],3]
    res = rbind(m_observed - expected_m, r_observed - expected_r)
    -0.5 * sum(rowSums((t(res) %*% ihalfSigma)^2))
  }, seq_len(length(all_data)))
  #log-sum-exp-trick
  unnorm_prob = exp(theta_prob_log - max(theta_prob_log))
  norm_prob = unnorm_prob / sum(unnorm_prob)
  sampled_theta = theta_grid[sample(1:nrow(grid_theta), 1, prob = norm_prob), ]
  theta_store[iter, ] = sampled_theta
  theta = sampled_theta
  
  # dr.pati suggestion:
  # x <- theta_prob_log * (theta_prob_log > -200000)
  # y <- exp(x - max(x))
  # res <- y/sum(y) ; print(res)
  
  # #gumbel-trick:
  # E = rexp(length(theta_prob_log), 1)  ; G = -log(E) # Generate Gumbel noise
  # gumbel_perturbed = theta_prob_log + G # Add Gumbel noise to the log probabilities
  # selected_index = which.max(gumbel_perturbed)  # Find the index of the maximum Gumbel-perturbed value
  # sampled_theta = theta_grid[selected_index, ]  # Sample the corresponding theta
  # theta_store[iter, ] = sampled_theta   # Store and update the theta
  # theta = sampled_theta
  
  #update sigma:
  #find residuals based on previously sampled theta and p
  theta_idx = which(apply(theta_grid, 1, function(row) all(row == sampled_theta)))
  row_indexx = match(p,all_data[[theta_idx]][,1])
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




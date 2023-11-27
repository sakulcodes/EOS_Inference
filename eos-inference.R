library(LaplacesDemon)
library(MASS) 
library(ggplot2)

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
correlation <- cor(observed_M, observed_R); print(correlation)
plot(observed_R, observed_M, xlab = "Observed Radius", ylab = "Observed Mass", main ="Observed Data")
obs_cov_mat <- cov(combined_data) ; obs_cov_mat#covariance matrix

#generate data:
n = length(observed_M)
correlated_noise = mvrnorm(n, mu = c(0, 0), Sigma = obs_cov_mat)
m_observed = observed_M + correlated_noise[,2] ; r_observed = observed_R + correlated_noise[,1]
plot(r_observed, m_observed, main ="Generated Data")

#Bayesian Parametric Model 
# (m_i, r_i)' | p_i, theta, Sigma  ~ N_2( {m(p_i, theta),r(p_i, theta)} , Sigma) 
#              theta  ~ Discrete Uniform(theta_1,....,theta_K) 
#               p_i | theta ~ Discrete Uniform {Pmin(theta),....,Pmax(theta)} independently for i = 1 to n
#               Sigma ~ Inverse_Wishart(nu = 2, Psi = diag(2))


# ------------------------ ORACLE MATRIX  - {M(pi,theta), R(pi,theta)} - Implicit Relationship - Solves the EOS & TOV:
#The FORTRAN oracle should solve (M(pi,theta), R(pi, theta)) FOR EVERY COMBINATION in the grid values for theta and p. 
#Imagine a N by K matrix, N = no. of central densities/central pressure, K = no.of theta 

#creating the oracle matrix:
num_files <- 100  # Assuming you have 100 files
base_path <- "/Users/sakul/Desktop/ldmcrust/res_eos_pnm32/eft_pnm32_"

# Step 1: Read Data and Create a Unique List of Central Densities
all_data <- vector("list", num_files)
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)[, 1:3]  # Read only the first three columns
  all_data[[i]] <- data_eos
}

# Step 2: Extract Unique Central Densities
unique_densities <- unique(unlist(lapply(all_data, function(x) x[,1])))
unique_densities <- sort(unique_densities, decreasing = TRUE)

# Step 3: Build the Oracle Matrices for Mass (M_Oracle) and Radius (R_Oracle)
M_Oracle <- R_Oracle <- matrix(NA, nrow = length(unique_densities), ncol = num_files)
rownames(M_Oracle) <- rownames(R_Oracle) <- unique_densities

for (i in 1:num_files) {
  for (j in 1:nrow(all_data[[i]])) {
    density <- all_data[[i]][j, 1]
    mass <- all_data[[i]][j, 2]
    radius <- all_data[[i]][j, 3]
    row_index <- which(rownames(M_Oracle) == density)
    M_Oracle[row_index, i] <- mass
    R_Oracle[row_index, i] <- radius
  }
}

M_Oracle ; dim(M_Oracle)
R_Oracle ; dim(R_Oracle)


# -----> Parametric curve visualization 
plot( R_Oracle[,1], M_Oracle[,1], xlab = "Radius", ylab = "Mass", col = "blue") ; points(R_Oracle[,2], M_Oracle[,2], col = "red")
points(R_Oracle[,3], M_Oracle[,3], col = "green")


#Check what linear combination results in the ({ai},{bi}) grid.
#grid for theta: #
prior <- read.table("/Users/sakul/Desktop/ldmcrust/data/eft_pnm32_kfr16_par_all.dat")[, -c(5, 6)]
grid_theta <- prior ; grid_theta <- grid_theta[1:100,]  
#for pseudo purposees - taking grid for only the first 100 values as only 100 data points are loaded in the oracle matrix. 



#THE FIRST COLUMN OF nsmr.dat file is the central densities

#grid for central density (nc)
all_densities <- c()
for (i in 1:num_files) {
  file_name <- sprintf("%s%06d_ldm_nsmr.dat", base_path, i)
  data_eos <- read.table(file_name)
  all_densities <- c(all_densities, data_eos$V1)
}
unique_densities <- sort(unique(all_densities), decreasing = TRUE) ; unique_densities
grid_p <- unique_densities ; log_prior_p = -log(length(grid_p))



# ------------------------------Gibbs-Sampling-Setup:
nreps = 500; nu = 2; Psi = diag(2) ; n = length(observed_M)
#storage variables:
theta_store = matrix(0, nrow = nreps, ncol = 8) ; p_store = matrix(0, nrow = nreps, ncol = n) ; Sigma_store = array(0, dim = c(2, 2, nreps))
#Initialize theta (take the row closest to the mean)
theta_grid <- as.matrix(grid_theta)
mean = colMeans(theta_grid) ; distances = apply(theta_grid, 1, function(row) sum((row - mean)^2))  # Calculate squared distances
theta = theta_grid[which.min(distances), ]

#Initialize Sigma:
Sigma = rinvwishart(nu, Psi)  #initialization

for (iter in 1:nreps) {
  iSigma = solve(Sigma)
  
  #Update p
  theta_idx = which(apply(theta_grid, 1, function(row) all(row == theta)))
  for (i in 1:n) {
    p_prob_log = mapply(function(p_val) {
      p_index = which(grid_p == p_val) 
      indicator = ifelse(is.na(M_Oracle[p_index, theta_idx]), 0, 1)  # Indicator function.
      if (indicator == 0) {
        return(0)
      } else {
        expected_m = M_Oracle[p_index, theta_idx]
        expected_r = R_Oracle[p_index, theta_idx]
        res = c(m_observed[i] - expected_m, r_observed[i] - expected_r)
        return(-0.5 * res %*% iSigma %*% res)
      }
    }, grid_p)
    unnorm_prob = exp(p_prob_log - max(p_prob_log))
    norm_prob = unnorm_prob / sum(unnorm_prob)
    p_store[iter, i] = sample(theta_grid, 1, prob = norm_prob)
  }
  p = p_store[iter,] #sampled pressure
  sampled_p_indexes = match(p, grid_p)
  
  #Update theta vector
  theta_prob_log = sapply(1:ncol(theta_grid), function(i) {  # Loop over columns (theta values)
    expected_m = M_Oracle[, i] ; expected_r = R_Oracle[, i]
    res = rbind(m_observed - expected_m[sampled_p_indexes], r_observed - expected_r[sampled_p_indexes])
    -0.5 * sum(diag(t(res) %*% iSigma %*% res))
  })
  unnorm_prob = exp(theta_prob_log - max(theta_prob_log))
  norm_prob = unnorm_prob / sum(unnorm_prob)
  sampled_row_idx = sample(1:nrow(theta_grid), 1, prob = norm_prob)
  sampled_theta = theta_grid[sampled_row_idx, ]
  theta_store[iter, ] = sampled_theta
  
  #update sigma: 
  #find residuals based on previously sampled theta and p
  theta_idx = which(apply(theta_grid, 1, function(row) all(row == sampled_theta)))
  p_index = which(grid_p == p_val)
  residuals_m = m_observed - M_Oracle[p_index, theta_idx] ; residuals_r = r_observed - R_Oracle[p_index, theta_idx] 
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



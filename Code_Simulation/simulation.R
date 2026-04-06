source("data_generator.R") # Generates X and Y =  X %*% beta + epsilon
source("loglikelihood.R")
source("noise_distribution.R")

library(sfsmisc) # to use function list_()
library(VGAM) # to use function rlaplace()

################################################################################
#                        Settings of the Simulation
################################################################################
# Settings of error distribution
noise <- "Gaussian"  # Choose between "Gaussian" and "Laplace"
sdError <- 0.8*sqrt(10) # standard deviation of noise


# Settings of true parameters
# Corners & centers of faces of the cube [-1,1]^3
pts <- rbind(c(2,2,2),
						 as.matrix(unname(expand.grid(rep(list(c(-1,1)),3)))),
						 c(1,0,0), c(-1,0,0), c(0,1,0),c(0,-1,0),c(0,0,-1), c(0,0,1))
# Small fluctuations to the points
set.seed(2025) 
pts_fl <- t(apply(pts, 1, function(x){x + rnorm(3,0,0.1)}))
nj <- nrow(pts) # number of SNRs
# Generate true parameters
beta0s <- pts_fl/apply(pts_fl,1,function(v){norm(v,"2")})*seq(0,8,length.out=nj)


# Settings of the covariate variable X
p <- ncol(beta0s) # dimension of X
X_distr <- "Gaussian" # Choose between "Gaussian" and "Uniform"
meanX <- rep(0,p) # mean vector of X
SigmaX <- diag(1,p) # covariance matrix of X


#X_distr <- "Gaussian" # Choose between "Gaussian" and "Uniform"
#meanX <- rep(5,p) # mean vector of X
#SigmaX <- diag(1,p) # covariance matrix of X

# X ~ U([-\sqrt{3},\sqrt{3}]^3)
#X_distr <- "Uniform" # Choose between "Gaussian" and "Uniform"
#meanX <- rep(0,p) # mean vector of X
#SigmaX <- diag(1,p) # covariance matrix of X

# X ~ U([5-\sqrt{3}, 5+\sqrt{3}]^3)
#X_distr <- "Uniform" # Choose between "Gaussian" and "Uniform"
#meanX <- rep(5,p) # mean vector of X
#SigmaX <- diag(1,p) # covariance matrix of X


# Settings of the simulation
ratios <- c(0.2, 0.6) # number of matched data / unmatched data
n_us <- c(200,1000,5000) # number of unmatched data
iter <- 2500 # number of iteration
save_results <- TRUE # Save the simulation results as RDS file? 
# TRUE: save as RDS file
# FALSE: do not save

################################################################################
#         Start of Simulation (No setting for the following codes)
################################################################################
# Check if the noise distribution is available
stopifnot((noise == "Gaussian") | (noise == "Laplace"))

for(rr in 1:length(ratios)){
	ratio <- ratios[rr]
	for (nn in 1:length(n_us)){
		n_u <- n_us[nn]
		n_m <- n_u * ratio # number of matched data
		
		
		##################################
		
		
		# results storage for all SNRs
		OLSEs <- vector("list", nj) # list to store OLSE results for each SNR
		MLEs_matched <- vector("list", nj) # MLE of matched samples
		MLEs <- vector("list", nj) # SSLEMLE
		
		# Simulation
		# Iteration for each SNR
		for (j in 1:nj){
			beta0 <- beta0s[j,] # true beta
			
			# results storage for current SNR
			OLSE <- matrix(rep(NA,iter*p), nrow = iter) # results of OLSE
			MLE_matched <- matrix(rep(NA,iter*p), nrow = iter) # results of matched MLE
			MLE <- matrix(rep(NA,iter*p), nrow = iter) # results of SSLEMLE
			
			for(i in 1:iter){
				
				set.seed(2025+i) # random seed
				
				# Generate matched data
				data_m <- data_generator(noise_type = noise, x_type = X_distr, 
																 dim = p, beta = beta0, num = n_m, 
																 mean = meanX, sigma = SigmaX) 
				x_m <- data_m[[1]]
				y_m <- data_m[[2]]
				
				# Generate unmatched covariates
				data_ux <- data_generator(noise_type = noise, x_type = X_distr,
																	dim = p, beta = beta0, num = n_u, 
																	mean = meanX, sigma = SigmaX) 
				# Generate unmatched responses
				data_uy <- data_generator(noise_type = noise, x_type = X_distr, 
																	dim = p, beta = beta0, num = n_u, 
																	mean = meanX, sigma = SigmaX) 
				x_u <- data_ux[[1]]
				y_u <- data_uy[[2]]
				
				################################### OLSE ###################################
				
				OLSesti <- unname(lm(y_m ~ x_m-1)$coefficients) # no intercept
				OLSE[i,] <- OLSesti
				
				################################### MLE ####################################
				# MLE of matched samples
				MLesti_matched <- optim(OLSesti,
																fn <- function(beta){
																	(-ll_match(noise, x_m, y_m, beta))},
																method = "BFGS")
				MLE_matched[i,] <- MLesti_matched$par
				
				# SSLEMLE
				MLesti <- optim(OLSesti,
												fn <- function(beta){
													(-ll_unmatch(noise, x_u, y_u, beta)
													 -ll_match(noise, x_m, y_m, beta))}, 
												method = "BFGS")
				MLE[i,] <- MLesti$par
				
			
				cat("Iteration",i,'/',iter, ';', 'parameter',  j, '/', nj, ';', 
						'unmatched data', nn ,'/', length(n_us), ';', 'ratios', rr, '/', length(ratios),
						'is done','\n')
			}
			
			OLSEs[[j]] <- OLSE
			MLEs_matched[[j]] <- MLE_matched
			MLEs[[j]] <- MLE
			
		}
		
		output <- list_ (OLSEs, MLEs_matched, MLEs, beta0s, 
										 sdError, noise, X_distr, meanX, SigmaX, iter,
										 n_m, n_u)
		
		if (save_results == TRUE){
			saveRDS(output, 
							file = paste0("simulation_dim=", p, 
														"_n_u=", n_u, "_ratio=", ratio,
														"_iter=", iter, "_noise=",noise, "_X=", X_distr,
														"_mean1=", meanX[1]))
		}
		
		
		
		##################################
		
		
		
	}
}


################################################################################
#            End of Simulation (No setting for the above codes)
################################################################################


library(readxl) # read Excel data
library(sfsmisc) # use list_()
library(VGAM) # use rlaplace()
source("loglikelihood.R") # log-likelihood functions
source("noise_distribution.R") # error distributions ("Gaussian" or "Laplace")

# Load data set
ccpp <- read_excel("Folds5x2_pp.xlsx") # Combined Cycle Power Plant
ccpp <- data.frame(ccpp)

################################################################################
#                           pre-process of data
################################################################################

# OLS based on whole data set (model without intercept after standardization)
fit <- lm(PE ~ ., 
					data = ccpp) # only consider AT and V
# Analysis of linear regression
summary(fit)
plot(fit,1) # TA-plot
plot(fit,2) # QQ-plot
# residuals 
hist(fit$residuals, probability = TRUE, nclass = 90, main = "Histogram of Residuals") 
lines(seq(-20,20,0.0001), dnorm(seq(-20,20,0.0001), 0, sigma(fit)),
			type = 'l', col = "blue")

################################################################################
#                          Setting of Simulation
################################################################################
n_u <- 1600 # number of unmatched data
ratio <- 10/n_u # m/n
n_m <- n_u * ratio # number of matched data
iter <- 1000 # number of iterations
noise <- "Gaussian" # we assume that the noise is Gaussian distributed
known_sdError <- 1 # use estimated s.d. of noise from OLS based on whole data (TRUE)
                   # OR use s.d. from OLS based on matched data (FALSE)

save_output <- 1 # save output as RDS file (TRUE)


seed <- 2025

################################################################################
#                   Start of Simulation (no more settings)
################################################################################
# Assume that the OLS based on all data points is the true parameter.
beta0 <- unname(fit$coefficients)
p <- length(beta0) # dimension of true beta

# results storage
OLSE <- matrix(rep(NA,iter*p), nrow = iter) # results of OLSE
MLE <- matrix(rep(NA,iter*p), nrow = iter) # results of pseudo MLE
mMLE <- matrix(rep(NA, iter*p), nrow = iter) # results of matched MLE
sigma_matched <- rep(NA, iter) # estimated s.d. of noise from OLS based on matched data

mse_OLSE <- rep(NA,iter)
mse_MLE <- rep(NA,iter)
mse_mMLE <- rep(NA,iter)


# Simulation
for (i in 1:iter){
	set.seed(seed+i)
	
	# divide whole data into testing set and training set
	
	index_test <- sample(nrow(ccpp), nrow(ccpp)/4) # index of testing set
	index_train <- setdiff(1:nrow(ccpp), index_test) # index of training set
	
	index_m <- sample(index_train, n_m) # sample matched data from training set
	index_ux <- sample(setdiff(index_train,index_m), n_u) # sample unmatched covariate data
	index_uy <- sample(setdiff(index_train,index_m), n_u) # sample unmatched response data
	
	# matched data based on random index
	x_m <- as.matrix(ccpp[index_m, -5])
	rownames(x_m) <- NULL
	x_m <- cbind(rep(1,n_m), x_m)
	y_m <- matrix(ccpp[index_m,5], ncol = 1)
	
	# unmatched data based on random index
	x_u <- as.matrix(ccpp[index_ux, -5])
	rownames(x_u) <- NULL
	x_u <- cbind(rep(1,n_u), x_u)
	y_u <- matrix(ccpp[index_uy,5], ncol = 1)
	
	# testing set
	x_test <- as.matrix(ccpp[index_test, -5])
	rownames(x_test) <- NULL
	x_test <- cbind(rep(1,length(index_test)),x_test)
	y_test <- matrix(ccpp[index_test,5], ncol = 1)

	################################### OLSE #####################################
	
	fit_matched <- lm(y_m ~ x_m-1) # OLS based on matched data
	sigma_matched[i] <- sigma(fit_matched)
	OLSesti <- unname(fit_matched$coefficients)
	OLSE[i,] <- OLSesti
	
	################  estimated standard deviation of noise  #####################

	sdError <- sigma(fit_matched) * (1-known_sdError) + sigma(fit) * known_sdError
	
	##############################  matched MLE ##################################	
	mMLesti <- optim(OLSesti,
									 fn <- function(beta){
									 	(-ll_match(noise, x_m, y_m, beta))}, 
									 method = "BFGS")$par
	
	mMLE[i,] <- mMLesti
	
	################################# SSLEMLE ####################################
	MLesti <- optim(OLSesti,
									fn <- function(beta){
										(-ll_unmatch(noise, x_u, y_u, beta)
										 -ll_match(noise, x_m, y_m, beta))}, 
									method = "BFGS")$par
	
	MLE[i,] <- MLesti
	
	#############################       MSE    ###################################
	
	mse_OLSE[i] <- mean((y_test - x_test %*% OLSesti)^2)
	mse_mMLE[i] <- mean((y_test - x_test %*% mMLesti)^2)
	mse_MLE[i] <- mean((y_test - x_test %*% MLesti)^2)
	
	##############################################################################
	
	cat("Iteration",i,'/',iter,'is done','\n')
	
}


output <- list_(OLSE, MLE, mMLE, mse_OLSE, mse_mMLE, mse_MLE,
								index_test, index_train, ccpp,
								n_m, n_u, iter, sigma_matched, known_sdError, noise)
if(save_output == TRUE){
	saveRDS(output, 
					file = paste0("ccpp_","n=",n_u,"_m=",n_m))
}
################################################################################
#                            End of Simulation
################################################################################
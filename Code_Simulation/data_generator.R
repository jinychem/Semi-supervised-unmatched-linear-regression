# Data generator of y = x %*% beta + epsilon
library(MASS)
data_generator <- function(noise_type, x_type, dim, beta, num, mean, sigma){
	p <- dim
	beta <- beta
	stopifnot((x_type == "Gaussian") | (x_type == "Uniform"))
	stopifnot(length(beta) == p)
	stopifnot(length(mean) == p)
	stopifnot((nrow(sigma) == p) && (ncol(sigma) == p))
	n <- num
	
	if (x_type == "Gaussian"){
		x <- mvrnorm(n = n, mu = mean, Sigma = sigma)
	}
	if (x_type == "Uniform"){
		stopifnot((p == 1) | (all(sigma == diag(diag(sigma)))))
		x <- matrix(rep(NA, n*p), ncol = p)
		for (k in 1:p){
			x[,k] <- runif(n, min = mean[k]- sqrt(3*diag(sigma)[k]), 
										 max = mean[k] + sqrt(3*diag(sigma)[k]))
		}
	}
	epsilon <- rdistr_epsilon(noise_type, n)
	y <- x %*% beta + epsilon
	return(list(x,y,epsilon))
}


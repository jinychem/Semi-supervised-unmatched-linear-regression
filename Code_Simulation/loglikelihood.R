ll_unmatch <- function(noise_type, x, y, beta) {
	a <- rep(NA,length(y))
	for (j in 1:length(y)){
		fy <- mean(ddistr_epsilon(noise_type, y[j] - x %*% beta))
		a[j] <- log(fy + 1e-16)
	}
	sum(a)
}



ll_match <- function(noise_type, x, y, beta){
	a <- rep(NA,length(y))
	for (k in 1:length(y)){
		a[k] <- log(ddistr_epsilon(noise_type, y[k] - x[k,] %*% beta) + 1e-16)
	}
	sum(a)
}
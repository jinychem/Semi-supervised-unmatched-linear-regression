# random sampling from given noise distribution
rdistr_epsilon <- function(noise_type, n){
	switch (noise_type,
					"Gaussian" = rnorm(n, mean = 0, sd = sdError), # epsilon ~ N(0,sdError^2)
					"Laplace" = rlaplace(n, location = 0, scale = sdError/sqrt(2))
	)
}
# CDF of noise distribution
pdistr_epsilon <- function(noise_type, n){
	switch(noise_type,
				 "Gaussian" = pnorm(n, mean = 0, sd = sdError), # epsilon ~ N(0,sdError^2)
				 "Laplace" = plaplace(n, location = 0, scale = sdError/sqrt(2))
	)
}
# density of noise distribution
ddistr_epsilon <- function(noise_type, x){
	switch (noise_type,
					"Gaussian" = dnorm(x, mean = 0, sd = sdError), # epsilon ~ N(0,sdError^2)
					"Laplace" = dlaplace(x, location = 0, scale = sdError/sqrt(2))
	)
}
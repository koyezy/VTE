############ Compute mean bias via toy simulation ############

library(statmod)

# use gauss.quad() from 'statmod' package
compute.bias.gh <- function(beta0, beta1, mu, sigma, c, n = 20) {
  
  # Gauss-Hermite nodes and weights
  gh <- gauss.quad(n, kind = "hermite")
  nodes <- gh$nodes # correspond to t
  weights <- gh$weights # correspond to w
  
  # Define expit
  expit = function(x) 1/(1+exp(-x))
  
  # Function g(t)
  g <- function(t) {
    y <- mu + sigma * sqrt(2) * t
    k <- beta0 + beta1*y
    {y*(1-expit(k)) + c*expit(k)}/ sqrt(pi)
  }
  
  # Gauss-Hermite approximation
  I_approx <- sum(weights * sapply(nodes, g)) %>% round(4)
  return(I_approx-mu)
}


find_k = function(input, target.rate){
  k = -1000
  expit = function(x){
    1/(1+exp(-x))
  }
  while(mean(expit(input+k)) <= target.rate){
    k = k + 0.01
  }
  return(k)
}

compute.bias <- function(mu, sigma, abs.rates, betas) {
  ns <- 10000
  results <- data.frame()
  
  for (abs.rate in abs.rates) {
    for (beta in betas) {
      # Generate outcome
      y <- rnorm(ns, mu, sigma)
      
      # Find intercept k for given beta and target absence rate
      k <- find_k(beta * y, abs.rate)
      
      # Missingness mechanism
      p <- 1 / (1 + exp(-(k + beta * y)))  # expit
      miss <- rbinom(ns, 1, p)
      
      # Single imputation with mean of observed values
      m <- median(y[miss == 0]) # c: imputed value
      ynew <- ifelse(miss == 1, m, y)
      
      # Compute mean bias
      bias <- mean(ynew - y)
      
      # Store results
      results <- rbind(
        results,
        data.frame(
          abs.rate = abs.rate,
          beta1 = beta,
          beta0 = round(k, 4),
          imputed.median = round(m, 4),
          bias.true = round(bias, 4),
          bias.mcar = round((m-mu)*abs.rate,4),
          bias.mar = compute.bias.gh(beta0=k, beta1=beta, mu, sigma, c=m, n = 20)
        )
      )
    }
  }
  
  return(results)
}

abs.rates <- c(0.1, 0.3)
# MCAR: no dependence between Y and M) when beta=0
# MAR: negative relationship between outcome and missing (pt with missing are less likley to develop y)
betas <- c(0, -0.2, -0.5)
res = compute.bias(mu = 0, sigma = 1, abs.rates = abs.rates, betas = betas)
print(res)

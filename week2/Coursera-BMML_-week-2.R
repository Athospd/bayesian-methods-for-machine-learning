dmvgaussian <- function(X, mu, sigma) {
  gaussian_normalizer = sqrt((2 * 3.14159265359)^length(mu) * det(sigma))
  
  x = t(X - mu)
  invSx = solve(sigma, x)
  xInvSx = colSums(-1/2 * (x * invSx))
  density = exp(xInvSx)/gaussian_normalizer
  
  return(density)
}

#' E_step
#'
#' @param X (N x d), data points
#' @param pi (C), mixture component weights 
#' @param mu (d x C), mixture component means
#' @param sigma (d x d x C), mixture component covariance matrices
#'
#' @return gamma: (N x C), probabilities of clusters for objects
#' 
#' @export
#'
#' @examples
E_step <- function(X, pi, mu, sigma) {
  
  N = nrow(X) # number of objects
  C = length(pi) # number of clusters
  d = ncol(mu) # dimension of each object
  gamma = matrix(0, N, C) # distribution q(T)
  
  for(c in seq(C))
    gamma[, c] = pi[c] * mvtnorm::dmvnorm(X, mu[c,], sigma[,,c])
  
  
  gamma = gamma/rowSums(gamma)
  
  return(gamma)
}

library(tidyverse)

samples <- read_rds("week2/samples.rds")

E_step(
  samples$data,
  samples$pi0,
  samples$mu0,
  samples$sigma0
)

b <- dmvgaussian(
  samples$data,
  samples$mu0[1,],
  samples$sigma0[,,1]
)

a <- mvtnorm::dmvnorm(
  samples$data,
  samples$mu0[1,],
  samples$sigma0[,,1]
)

hist(a - b)

library(ggplot2)
library(dplyr)

data.frame(
  x = rexp(1000) 
) %>%
  ggplot(aes(x = x)) + geom_density(colour = "red")

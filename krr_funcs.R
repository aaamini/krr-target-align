# Needs "p" to be defined
h <- sqrt(p/2) #.1

Kfun <- function(x, y, sig2=h^2)   exp(-sum((x-y)^2)/(2*sig2))

create_kernel_matrix = function(X,Y) {
  apply(Y, 1, function(y) apply(X,  1, function(x) Kfun(x,y)))
  # 1 + apply(Y, 1, function(y) apply(X,  1, function(x) Kfun(x,y)))
}

create_true_fun = function(n_knots = 10, p, coeffs = rnorm(n_knots)) {
  X_knots = matrix(runif(n_knots*p), nrow = n_knots)
  
  # x has to be a p x 1 vector
  
  function(X) create_kernel_matrix(X, X_knots) %*% coeffs
}

# *** Now implemented in C++
# get_theo_mse = function(lambda, v, mu_r, sig2_over_n) {
#   r = length(mu_r)
#   n = length(v)
#   Gamma_diag = rep(0, n)
#   Gamma_diag[1:r] = mu_r / (mu_r + lambda)
#   Gamma = Diagonal(x = Gamma_diag)
#   eye = Diagonal(n)
#   ex_approx_err = sum(((eye - Gamma) %*% v)^2)
#   ex_estim_err = sig2_over_n * sum(diag(Gamma^2))
#   ex_approx_err + ex_estim_err
# }

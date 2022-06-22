library(Matrix)
library(ggplot2)
library(dplyr)

set.seed(1234)
p <- 4  # input dimension
source("krr_funcs.R")

# The true function -- a linear combination of "n_knots" kernels functions
f_tru = create_true_fun(n_knots = 20, p = p)

n = 200
X = matrix(runif(n*p), nrow = n, ncol = p)

fs = f_tru(X)
sd_signal <- sd(fs)
sig = .75*sd_signal
y <- fs + rnorm(n, sd = sig)

r = 20
K = create_kernel_matrix(X,X) / n 
out = eigen(K)
U = out$vectors
U_r = U[ , 1:r]
mu_r = out$values[1:r]
D_r = diag(mu_r)
Kt = U_r %*% D_r %*% t(U_r) 

# compute empirical MSE
lambda_vec = 10^seq(-4,-2, length.out = 20)
runs = expand.grid(lambda = lambda_vec, rep = 1:50)
res = do.call(rbind, lapply(1:nrow(runs), function(j) {
  lambda = runs[j, "lambda"]
  rep = runs[j, "rep"]
  y <- fs + rnorm(n, sd = sig)
  ft = solve(Kt + lambda*diag(n), Kt %*% y)
  data.frame(lambda = lambda, n_mse = mean((ft - fs)^2), rep = rep)
}))

res2 = res %>% 
  group_by(lambda) %>% 
  summarise(n_mse = mean(n_mse)) %>% 
  mutate(type = "Empirical")

# compute the theoretical MSE
v = t(U) %*% fs / sqrt(n)

ex_n_mse = sapply(lambda_vec, function(lambda) {
  Gamma_diag = rep(0, n)
  Gamma_diag[1:r] = mu_r / (mu_r + lambda)
  Gamma = Diagonal(x = Gamma_diag)
  eye = Diagonal(n)
  ex_approx_err = sum(((eye - Gamma) %*% v)^2)
  ex_estim_err = (sig^2/n) * sum(diag(Gamma^2))
  ex_n_mse = ex_approx_err + ex_estim_err
  ex_n_mse  
})

res3 = bind_rows(res2, 
                 data.frame(lambda = lambda_vec, n_mse = ex_n_mse, type = "Theory"))

res3 %>% ggplot(aes(x = lambda, y = n_mse, color = type)) + 
  geom_line(size = 1.1) + 
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("Normalized MSE") + xlab("Lambda") 
# ggsave("mse_krr.png", width = 4, height = 4)

# components above the noise level
which(v^2 > sig^2/n)


     
rm(list=ls())
#
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(plotly)
library(reticulate)

folder = "."

ppathFigures <- file.path(folder,"simulation")
if (!file.exists(ppathFigures)) {
  dir.create(ppathFigures, recursive = T)
} 

p = 4
set.seed(1234)
source(file.path(folder,"krr_funcs.R"))
Rcpp::sourceCpp(file.path(folder,"src","krr_calc.cpp"), verbose = T)


# The true function -- a linear combination of "n_knots" kernels functions
# This is not used in this code ... just there to advance the RNG seed to match the plots in paper
f_tru = create_true_fun(n_knots = 20, p = p)

# Create Matrix X
n = 200
X = matrix(runif(n*p), nrow = n, ncol = p)

# Calculate the kernel matrix
K = create_kernel_matrix(X,X) / n 
# Spectrum of K
out = eigen(K)
U = out$vectors
mu = out$values

#l: where the alignment starts
#k: how many components are aligned
#r: where is the truncation
b=30
l=0

# Generate the signal using the alignment pattern VV
VV <- rep(0,n)
VV[1:b+l] <- rnorm(b, 0 , 1)
VV <- VV / sqrt(sum(VV^2))

#Let v==VV
fs <- sqrt(n) * U%*%VV
sd_signal <- sd(fs)
v = matrix(VV, nrow=length(VV), ncol=1)
#-----------------------------------------

#store results in res
res <- NULL
#different values of lambda, the regularization parameter of the kernel ridge regression
rseq <- seq(1, 15, by=1)
lambda_vec = 10^seq(-5, 1, length=5000)
sig = 2

# Computes the entire MSE surface in C++
ex_n_mse = get_theo_mse_3d(lambda_vec, rseq, sig^2/n, v, mu)

# # This is inefficient ... TODO: Create Dataframe directly in C++
# for (ri in seq_along(rseq)) {
#   for (li in seq_along(lambda_vec)) {
#     res = bind_rows(res, data.frame(lambda = lambda_vec[li],
#                                     n_mse = ex_n_mse[li, ri, 1],
#                                     sig = sig,
#                                     r = rseq[ri]))
#   }
# }

m <- list(
  l = 50,
  r = 50,
  b = 100,
  t = 70,
  pad = 4
)

comb <- paste0("l=", l, ", b=", b, ", sig=", sig)

fig <- plot_ly(x=-log(lambda_vec), y=rseq, z=t(ex_n_mse[,,1]), 
               type="contour", line = list(width = 0), ncontours = 25) %>% 
  layout(title = list(text=paste0("MSE for different lambda and sigma values \n", comb), 
                      font =list(size=20, face="bold", color='black')),
         plot_bgcolor = "#e5ecf6",
         xaxis = list(
           title=list(text="-log(lambda)"),
           titlefont = list(size = 15, color="black"), 
           tickfont=list(size=10),
           automargin=TRUE, standoff = 5), 
         yaxis = list(
           title=list(text='r'), 
           titlefont = list(size = 15, color="black"),
           tickfont = list(size=10),
           automargin=TRUE, standoff = 50), 
        legend = list(orientation = "h",   
                       xanchor = "center",  
                       y = -0.1),
         margin = m) %>% config(mathjax = 'cdn')

fig

save_image(fig, file.path(ppathFigures, "mse_rank_vs_lambda_contour.pdf"), width = 500)

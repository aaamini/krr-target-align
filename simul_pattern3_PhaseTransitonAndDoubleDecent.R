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
k=30
l=0
r=30

#get "r" eigen vectors and values
U_r = U[ , 1:r]
mu_r = mu[1:r]


# Generate the signal using the alignment pattern VV
VV <- rep(0,n)
VV[1:k+l] <- rnorm(k, 0 , 1)
VV <- VV / sqrt(sum(VV^2))

#Let v==VV
fs <- sqrt(n) * U%*%VV
sd_signal <- sd(fs)
v = matrix(VV, nrow=length(VV), ncol=1)
#-----------------------------------------


#different values of lambda, the regularization parameter of the kernel ridge regression
lambda_vec = 10^seq(-5, 1, length=4000)

#different noise levels
sig_vec <- seq(0, 4, by=0.2) * sd_signal


#store results in res
res <- NULL

    
for(sig in sig_vec){
    ex_n_mse = sapply(lambda_vec, function(lam) get_theo_mse(lam, v, mu_r, sig^2/n))
  
    res = bind_rows(res, 
                    data.frame(lambda = lambda_vec, n_mse = ex_n_mse, sig=sig, 
                               r=r))
    
}

res$sig2 = round(res$sig, 2)
res$lambda2 <- -log(res$lambda)

#combination of three parameters: k, l and r
comb <- paste0("l=", l, ", b=", k, ", r=", r)


#-----plot MSE vs. -log(lambda) separate noise levels------------- 
res2 <- res %>% filter(sig2 %in% c(0, 2.4, 2.8, 3, 3.4, 4))

plt <- res2 %>% ggplot(aes(x = lambda2, y = n_mse)) + 
  facet_wrap(~sig2, scale='free', labeller = label_bquote(sigma == .(sig2))) +  
  geom_line(size = 1.2) + 
  theme_minimal() + 
  ggplot2::theme(legend.background = ggplot2::element_blank(),
                 #legend.title = ggplot2::element_blank(),
                 legend.position = c(0.9, 0.2)) + 
  labs(col='Value of r') + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("MSE") + xlab(expression(-log(lambda))) + 
  #ggtitle(expression(paste("MSE vs. ", -log(lambda), " \n k=",k, " l=", l, " r=", r))) + 
  ggtitle(bquote(atop("MSE vs." -log(lambda), .(comb)))) + 
  theme(plot.title = element_text(size=20, hjust = 0.5))+
  theme(strip.text.x = element_text(size = 20, angle = 0))+
  theme(axis.text.x = element_text( size = 12 )) +
  theme(axis.text.y = element_text( size = 12 )) +
  theme(axis.title = element_text( size = 16 ))
  
ggsave(file.path(ppathFigures, "mseVSlambda_pattern3_panel.pdf"))

#---------------------------contour plot---------------
mses <- matrix(res$n_mse, nrow=length(lambda_vec), ncol=length(sig_vec))


m <- list(
  l = 50,
  r = 50,
  b = 100,
  t = 70,
  pad = 4
)

fig <- plot_ly(x=-log(lambda_vec), y=sig_vec, z=t(mses), 
               type="contour") %>% 
  layout(title = list(text=paste0("MSE for different lambda and sigma values \n", comb), 
                      font =list(size=20, face="bold", color='black')),
         plot_bgcolor = "#e5ecf6",
         xaxis = list(
           title=list(text="-log(lambda)"),
           titlefont = list(size = 15, color="black"), 
           tickfont=list(size=10),
           automargin=TRUE, standoff = 5), 
         yaxis = list(
           title=list(text='sigma'), 
           titlefont = list(size = 15, color="black"),
           tickfont = list(size=10),
           automargin=TRUE, standoff = 50), 
        legend = list(orientation = "h",   
                       xanchor = "center",  
                       y = -0.1),
         margin = m) %>% config(mathjax = 'cdn')

fig

save_image(fig, file.path(ppathFigures, "mseVSlambda_pattern3_contour.pdf"), width = 500)

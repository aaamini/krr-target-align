#==================================================================================
#                     There are two segments of alignments   #                    #
#                                                                                 #
#==================================================================================

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



#different values of lambda, the regularization parameter of the kernel ridge regression
lambda_vec = 10^seq(-5, 1, length=100) # length used to be 4000

#sample size
n=200
#n=600

#alignment width of each segment
k=10


X = matrix(runif(n*p), nrow = n, ncol = p)

# Calculate the kernel matrix
K = create_kernel_matrix(X,X) / n 
# Spectrum of K
out = eigen(K)
U = out$vectors
mu = out$values


#l: where the alignment starts
l1=0
l2=30

# Generate the signal using the alignment pattern VV
VV <- rep(0,n)
VV[1:k+l1] <- 1/k
VV[1:k+l2] <- 1/k
VV <- VV / sqrt(sum(VV^2))

#Let v==VV
fs <- sqrt(n) * U%*%VV
sd_signal <- sd(fs)
v = matrix(VV, nrow=length(VV), ncol=1)
#-------------------------------------------------
#TAKE SIGMA AS RATIO OF NOISE LEVEL VS SAMPLE SIZE 
#different noise levels
sig_vec <- seq(0, 0.3, by=0.01) *  sqrt(n)
#-------------------------------------------------



#store results in res
res <- NULL

#r: where is the truncation

rseq <- seq(1, 50, by=1)

# Computes the entire MSE surface in C++
ex_n_mse = get_theo_mse_3d(lambda_vec, rseq, sig_vec^2/n, v, mu)

# This is inefficient ... TODO: Create Data fame directly in C++
for (ri in seq_along(rseq)) {
  for (si in seq_along(sig_vec)) {
    res = bind_rows(res, data.frame(lambda = lambda_vec,
                                    n_mse = ex_n_mse[ , ri, si],
                                    sig = sig_vec[si],
                                    r = rseq[ri]))
  }
}

# for(r in rseq){
#   #get "r" eigen vectors and values
#   U_r = U[ , 1:r]
#   mu_r = mu[1:r]
# 
#   for(sig in sig_vec){
#     ex_n_mse = sapply(lambda_vec, function(lam) get_theo_mse(lam, v, mu_r, sig^2/n))
# 
#     res = bind_rows(res,
#                     data.frame(lambda = lambda_vec, n_mse = ex_n_mse, sig=sig,
#                                r=r))
# 
#   }
# 
# }

# res$sig2 <- as.factor(paste0("noise sd ratio=", res$sig))
res$sig2 = round(res$sig / sqrt(n), 3)
res$sig3 = paste0("sigma/sqrt(n)=", res$sig2)
res$lambda2 <- -log(res$lambda)
res$r2 <- factor(paste0("r=", res$r), levels=paste0("r=", rev(rseq)))

# save(res, file=paste0(ppathResults, 
#                       "/mseVSr-twoPieceAlign", 
#                       "_k",k,"_l1",l1, "_l2",l2,"_n", n, "_ratioSigmaToN.rda"))

#-------------------------------------------------------------------------------
#-------------------------------------Draw plots--------------------------------
# 
# rm(list=ls())
# library(dplyr)
# library(ggplot2)
# library(plotly)
# 
# folder <- 'D:/Abbvie2020/papers/Kernel_target alignment_Arash'
# ppathResults <- file.path(folder,"simulation","results/twoSegment/")


# n=200
# #n=600
# k=10; l1=0; l2=30
# load(file=paste0(ppathResults, "/mseVSr-twoPieceAlign", "_k",k,"_l1",l1, "_l2",l2,"_n", n, "_ratioSigmaToN.rda"))


lambda.pick <- 1
lambda.pick.val <- unique(res$lambda2)[lambda.pick]
res.pick <- res %>% filter(lambda2==lambda.pick.val)
    
comb <- paste0("k=", k, ", l1=", l1, ", l2=", l2, ", n=", n, ", -log(\\lambda)=", round(lambda.pick.val,2))
lambda <- round(lambda.pick.val,1)
#----------plot MSE vs. -log(lambda) overlay all noise levels together
plt <- res.pick %>% filter(sig2 %in% c(0, 0.06, 0.12, 0.18, 0.24, 0.3)) %>%
  ggplot(aes(x = r, y = n_mse)) + 
  facet_wrap(~sig2, scales="free", labeller = label_bquote(sigma/sqrt(n) == .(sig2))) + 
  geom_line(size = 1.2) + 
  theme_minimal() + 
  ggplot2::theme(legend.background = ggplot2::element_blank(),
                 #legend.title = ggplot2::element_blank(),
                 legend.position = "bottom") + 
  scale_colour_gradientn(name=expression(sigma), 
                         colours=scales::muted(rainbow(10), l = 70, c = 80)) +
  ylab("MSE") + xlab("r") + 
  theme(plot.title = element_text(size=20, hjust = 0.5))+
  theme(strip.text.x = element_text(size = 20, angle = 0))+
  theme(axis.text.x = element_text( size = 12 )) +
  theme(axis.text.y = element_text( size = 12 )) +
  theme(axis.title = element_text( size = 16 ))+ 
  ggtitle(bquote(atop("MSE vs. r", 
                      l[1]==.(l1)*","~l[2]==.(l2)*","~b==.(k)*","~
                      -log(lambda)==.(lambda))))
  
plt
    
#-------------------------------contour plot-------------------------
library(plotly)
mses <- matrix(res.pick$n_mse, 
               nrow=length(unique(res.pick$r)), 
               ncol=length(unique(res.pick$sig2)))

m <- list(
  l = 50,
  r = 50,
  b = 50,
  t = 100,
  pad = 4
)

comb <- 'l1=0, l2=30, b=10, -log(lambda)=11.5'

fig <- plot_ly(x=unique(res.pick$r), 
               y=unique(res.pick$sig2), 
               z=t(mses), 
               type="contour",
               # autocontour = T,
               line = list(width = 0)
               ) %>% 
  layout(title = list(text=paste0('MSE  for  different r and sigma values \n',
                      'l1=0, l2=30, b=10, -log(lambda)=11.5'), 
                      font =list(size=20, face="bold", color='black')),
         plot_bgcolor = "#e5ecf6",
         xaxis = list(
           title=list(text="r"),
           titlefont = list(size = 15, color="black"), 
           tickfont=list(size=10),
           automargin=TRUE, standoff = 15), 
         yaxis = list(
           title=list(text='sigma/sqrt(n)'), 
           titlefont = list(size = 15, color="black"),
           tickfont = list(size=10),
           automargin=TRUE, standoff = 50), 
         legend = list(orientation = "h",   
                       xanchor = "center",  
                       y = -0.1),
         margin = m) %>% config(mathjax = 'cdn')



save_image(fig, file.path(ppathFigures, "mseVSr_twoseg_contour.pdf"), width = 500)
#

    
    
    
    
    
   
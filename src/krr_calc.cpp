// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
// using namespace arma;   // causes Reference to RcppArmadillo is ambiguous

// [[Rcpp::export]]
arma::vec get_theo_mse(
    arma::vec& lambda, 
    arma::vec& v,  // Full sequence of target alignment scores 
    arma::vec& mu,  // Truncated sequence of eigenvalues
    double sig2_over_n)
{
    int r = mu.n_elem;
    int n = v.n_elem;
    int n_lam = lambda.n_elem;
    arma::vec mse(n_lam, arma::fill::zeros);

    for (int l = 0; l < n_lam; l++) 
    {
      mse(l) = 0;
      for (int i = 0; i < r; i++) 
      {
        double gam_i = mu(i) / (mu(i) + lambda(l));
        mse(l) += pow((1-gam_i)*v(i), 2) + pow(gam_i, 2)*sig2_over_n;
      }
      for (int i = r; i < n; i++) {
        mse(l) += pow(v(i), 2);
      }
    } 
    return mse;
} 

// [[Rcpp::export]]
arma::cube get_theo_mse_3d(
    arma::vec& lambda, 
    arma::vec& rvec,
    arma::vec& sig2_over_n,
    arma::vec& v,  // Full sequence of target alignment scores 
    arma::vec& mu  // Full sequence of eigenvalues
    )
{
  int n_r = rvec.n_elem;
  int n_lam = lambda.n_elem;
  int n_sig = sig2_over_n.n_elem;
  arma::cube mse(n_lam, n_r, n_sig, arma::fill::zeros);
  // Rcpp::Rcout << mse;

  for (int si = 0; si < n_sig; si++) {
    for (int ri = 0; ri < n_r; ri++) {
      // Rcpp::Rcout << mu.rows(0, rvec(ri)-1);
      arma::vec mu_truncated = mu.rows(0, rvec(ri)-1);
      mse.slice(si).col(ri) = get_theo_mse(lambda, v, mu_truncated, sig2_over_n(si));
    }
  }
  return mse;
}
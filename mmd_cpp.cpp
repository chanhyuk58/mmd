#include <Rcpp.h>
using namespace Rcpp;

// Full Sample Objective Function Q
// [[Rcpp::export]]
double Q_obj_cpp(const NumericMatrix & x_mat,
                 const NumericVector & v0_vec,
                 const NumericVector & v1_vec,
                 const NumericVector & params,
                 const NumericVector & eta_vec) {

  int n = x_mat.nrow();
  int d = x_mat.ncol();
  
  // Use raw pointers for maximum speed during optimization loops
  const double *xptr = x_mat.begin();
  const double *v0ptr = v0_vec.begin();
  const double *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin();
  const double *pptr = params.begin();

  double gamma0 = pptr[0];
  double gammaV = pptr[1];
  double sumQ = 0.0;

  for (int i = 0; i < n; ++i) {
    double fx = 0.0;
    for (int j = 0; j < d; ++j) {
      // x_mat is column-major in memory: index is i + j * n
      fx += xptr[i + j * n] * pptr[2 + j];
    }
    
    double f1 = gamma0 + gammaV * v1ptr[i] + fx;
    double f0 = gamma0 + gammaV * v0ptr[i] + fx;
    double eta = etaptr[i];

    // The "kinks" (non-differentiable points) of the MMD objective
    if (f1 < eta) sumQ += (eta - f1);
    if (f0 > eta) sumQ += (f0 - eta);
  }

  return sumQ / static_cast<double>(n);
}

// Subsample Objective Function (Following Kaido, Molinari, and Stoyse 2019)
// [[Rcpp::export]]
double Q_obj_subsample_cpp(const NumericMatrix & x_mat,
                           const NumericVector & v0_vec,
                           const NumericVector & v1_vec,
                           const NumericVector & params,
                           const NumericVector & eta_vec,
                           const IntegerVector & idx_vec) {

  int b = idx_vec.size(); // Subsample size
  int n = x_mat.nrow();   // Full sample size (needed for matrix indexing)
  int d = x_mat.ncol();
  
  const double *xptr = x_mat.begin();
  const double *v0ptr = v0_vec.begin();
  const double *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin();
  const double *pptr = params.begin();
  const int *idxptr = idx_vec.begin();

  double gamma0 = pptr[0];
  double gammaV = pptr[1];
  double sumQ = 0.0;

  for (int k = 0; k < b; ++k) {
    int i = idxptr[k]; // Get the actual row index for this subsample draw
    
    double fx = 0.0;
    for (int j = 0; j < d; ++j) {
      fx += xptr[i + j * n] * pptr[2 + j];
    }
    
    double f1 = gamma0 + gammaV * v1ptr[i] + fx;
    double f0 = gamma0 + gammaV * v0ptr[i] + fx;
    double eta = etaptr[i];

    if (f1 < eta) sumQ += (eta - f1);
    if (f0 > eta) sumQ += (f0 - eta);
  }

  return sumQ / static_cast<double>(b);
}

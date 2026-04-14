#include <Rcpp.h>

using namespace Rcpp;

// Full Sample Objective Function
// [[Rcpp::export]]
double Q_obj_cpp(const NumericMatrix & x_mat,
                 const NumericVector & v0_vec,
                 const NumericVector & v1_vec,
                 const NumericVector & params,
                 const NumericVector & eta_vec) {

  int n = x_mat.nrow();
  int d = x_mat.ncol();
  const double *xptr = x_mat.begin();
  const double *v0ptr = v0_vec.begin();
  const double *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin();
  const double *pptr = params.begin();

  double g0 = pptr[0];
  double gV = pptr[1];
  double sumQ = 0.0;

  for (int i = 0; i < n; ++i) {
    double fx = 0.0;
    for (int j = 0; j < d; ++j) {
      fx += xptr[i + j * n] * pptr[2 + j];
    }
    
    double f1 = g0 + gV * v1ptr[i] + fx;
    double f0 = g0 + gV * v0ptr[i] + fx;
    
    double low = (f1 < f0) ? f1 : f0;
    double high = (f1 < f0) ? f0 : f1;
    double eta = etaptr[i];

    if (eta < low) {
      double diff = low - eta;
      sumQ += diff * diff;
    } else if (eta > high) {
      double diff = eta - high;
      sumQ += diff * diff;
    }
  }
  
  return sumQ / static_cast<double>(n);
}

// Subsample Objective Function
// [[Rcpp::export]]
double Q_obj_subsample_cpp(const NumericMatrix & x_mat,
                           const NumericVector & v0_vec,
                           const NumericVector & v1_vec,
                           const NumericVector & params,
                           const NumericVector & eta_vec,
                           const IntegerVector & idx_vec) {

  int b = idx_vec.size();
  int n = x_mat.nrow();
  int d = x_mat.ncol();
  
  const double *xptr = x_mat.begin();
  const double *v0ptr = v0_vec.begin();
  const double *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin();
  const double *pptr = params.begin();
  const int *idxptr = idx_vec.begin();

  double g0 = pptr[0];
  double gV = pptr[1];
  double sumQ = 0.0;

  for (int k = 0; k < b; ++k) {
    int i = idxptr[k]; 
    
    double fx = 0.0;
    for (int j = 0; j < d; ++j) {
      fx += xptr[i + j * n] * pptr[2 + j];
    }
    
    double f1 = g0 + gV * v1ptr[i] + fx;
    double f0 = g0 + gV * v0ptr[i] + fx;
    
    double low = (f1 < f0) ? f1 : f0;
    double high = (f1 < f0) ? f0 : f1;
    double eta = etaptr[i];

    if (eta < low) {
      double diff = low - eta;
      sumQ += diff * diff;
    } else if (eta > high) {
      double diff = eta - high;
      sumQ += diff * diff;
    }
  }

  return sumQ / static_cast<double>(b);
}

#include <Rcpp.h>
using namespace Rcpp;

// Full Sample Objective Function (Q)
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

  double g0 = pptr[0], gV = pptr[1], sumQ = 0.0;

  for (int i = 0; i < n; ++i) {
    double fx = 0.0;
    for (int j = 0; j < d; ++j) fx += xptr[i + j * n] * pptr[2 + j];
    
    double f1 = g0 + gV * v1ptr[i] + fx;
    double f0 = g0 + gV * v0ptr[i] + fx;
    
    double low = (f1 < f0) ? f1 : f0;
    double high = (f1 < f0) ? f0 : f1;
    double eta = etaptr[i];

    if (eta < low) {
      double d = low - eta; sumQ += d * d;
    } else if (eta > high) {
      double d = eta - high; sumQ += d * d;
    }
  }
  return sumQ / static_cast<double>(n);
}

// Subsample Objective Function (following Kaido, Molinari, and Stoye 2019)
// [[Rcpp::export]]
double Q_obj_subsample_cpp(const NumericMatrix & x_mat,
                           const NumericVector & v0_vec,
                           const NumericVector & v1_vec,
                           const NumericVector & params,
                           const NumericVector & eta_vec,
                           const IntegerVector & idx_vec) {

  int b = idx_vec.size(), n = x_mat.nrow(), d = x_mat.ncol();
  const double *xptr = x_mat.begin(), *v0ptr = v0_vec.begin(), *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin(), *pptr = params.begin();
  const int *idxptr = idx_vec.begin();

  double g0 = pptr[0], gV = pptr[1], sumQ = 0.0;

  for (int k = 0; k < b; ++k) {
    int i = idxptr[k];
    double fx = 0.0;
    for (int j = 0; j < d; ++j) fx += xptr[i + j * n] * pptr[2 + j];
    
    double f1 = g0 + gV * v1ptr[i] + fx;
    double f0 = g0 + gV * v0ptr[i] + fx;
    double low = (f1 < f0) ? f1 : f0, high = (f1 < f0) ? f0 : f1, eta = etaptr[i];

    if (eta < low) {
      double d = low - eta; sumQ += d * d;
    } else if (eta > high) {
      double d = eta - high; sumQ += d * d;
    }
  }
  return sumQ / static_cast<double>(b);
}

// Grid Search Evaluator
// [[Rcpp::export]]
NumericVector Q_eval_grid_cpp(const NumericMatrix & x_mat,
                              const NumericVector & v0_vec,
                              const NumericVector & v1_vec,
                              const NumericMatrix & param_grid,
                              const NumericVector & eta_vec) {
  
  int n_grid = param_grid.nrow(), p = param_grid.ncol();
  int n = x_mat.nrow(), d = x_mat.ncol();
  NumericVector results(n_grid);
  
  const double *xptr = x_mat.begin(), *v0ptr = v0_vec.begin(), *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin(), *gridptr = param_grid.begin();

  for (int g = 0; g < n_grid; ++g) {
    // Extract params for this grid point (R matrices are col-major)
    double g0 = gridptr[g + 0 * n_grid];
    double gV = gridptr[g + 1 * n_grid];
    
    double sumQ = 0.0;
    for (int i = 0; i < n; ++i) {
      double fx = 0.0;
      for (int j = 0; j < d; ++j) fx += xptr[i + j * n] * gridptr[g + (2 + j) * n_grid];
      
      double f1 = g0 + gV * v1ptr[i] + fx;
      double f0 = g0 + gV * v0ptr[i] + fx;
      double low = (f1 < f0) ? f1 : f0, high = (f1 < f0) ? f0 : f1, eta = etaptr[i];

      if (eta < low) {
        double d = low - eta; sumQ += d * d;
      } else if (eta > high) {
        double d = eta - high; sumQ += d * d;
      }
    }
    results[g] = sumQ / static_cast<double>(n);
  }
  return results;
}

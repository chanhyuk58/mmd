// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// General Linear MMD Logic
inline double compute_mmd_linear(const double* xptr, 
                                 const double* v0ptr, 
                                 const double* v1ptr, 
                                 const double* etaptr,
                                 const double* gridptr, 
                                 int g, int n_grid, 
                                 int n, int d) {
  
  double g0 = gridptr[g + 0 * n_grid];
  double gV = gridptr[g + 1 * n_grid];
  double sumQ = 0.0;

  for (int i = 0; i < n; ++i) {
    double fx = 0.0;
    for (int j = 0; j < d; ++j) {
      fx += xptr[i + j * n] * gridptr[g + (2 + j) * n_grid];
    }
    
    double f1 = g0 + gV * v1ptr[i] + fx;
    double f0 = g0 + gV * v0ptr[i] + fx;
    
    double low = (f1 < f0) ? f1 : f0;
    double high = (f1 < f0) ? f0 : f1;
    double eta = etaptr[i];

    // Squared Distance Penalty
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

// ==========================================================================
// FUNCTIONS
// ==========================================================================

// [[Rcpp::export]]
double Q_obj_cpp(const NumericMatrix & x_mat, const NumericVector & v0_vec,
                 const NumericVector & v1_vec, const NumericVector & params,
                 const NumericVector & eta_vec) {
  return compute_mmd_linear(x_mat.begin(), v0_vec.begin(), v1_vec.begin(), 
                            eta_vec.begin(), params.begin(), 0, 1, x_mat.nrow(), x_mat.ncol());
}

// [[Rcpp::export]]
// Projection method proposed by Kaido, Molinari, and Stoye 2019
double Q_obj_subsample_cpp(const NumericMatrix & x_mat, const NumericVector & v0_vec,
                           const NumericVector & v1_vec, const NumericVector & params,
                           const NumericVector & eta_vec, const IntegerVector & idx_vec) {
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
      double diff = low - eta; sumQ += diff * diff;
    } else if (eta > high) {
      double diff = eta - high; sumQ += diff * diff;
    }
  }
  return sumQ / static_cast<double>(b);
}

// [[Rcpp::export]]
NumericVector Q_eval_grid_cpp(const NumericMatrix & x_mat, const NumericVector & v0_vec,
                              const NumericVector & v1_vec, const NumericMatrix & param_grid,
                              const NumericVector & eta_vec) {
  int n_grid = param_grid.nrow(), n = x_mat.nrow(), d = x_mat.ncol();
  NumericVector results(n_grid);
  
  const double *xptr = x_mat.begin(), *v0ptr = v0_vec.begin(), *v1ptr = v1_vec.begin();
  const double *etaptr = eta_vec.begin(), *gridptr = param_grid.begin();

  // OMP Implementation if possible
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int g = 0; g < n_grid; ++g) {
    results[g] = compute_mmd_linear(xptr, v0ptr, v1ptr, etaptr, gridptr, g, n_grid, n, d);
  }
  return results;
}

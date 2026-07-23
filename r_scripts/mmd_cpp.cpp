#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

enum ModelFamily {
  FAMILY_LPM    = 0,
  FAMILY_LOGIT  = 1
};

inline double inv_logit(const double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  } else {
    const double z = std::exp(x);
    return z / (1.0 + z);
  }
}

inline double mean_function(const double lp, const ModelFamily family) {
  switch (family) {
  case FAMILY_LPM:
    return lp;
  case FAMILY_LOGIT:
    return inv_logit(lp);
  default:
    Rcpp::stop("Unknown model family.");
  }
  return NA_REAL;
}

double compute_objective(const NumericMatrix &x_mat,
                         const NumericVector &v0_vec,
                         const NumericVector &v1_vec,
                         const NumericVector &params,
                         const NumericVector &eta_vec,
                         const IntegerVector *idx,
                         const ModelFamily family) {

  const int n = (idx == nullptr) ? x_mat.nrow() : idx->size();
  const int N = x_mat.nrow();
  const int d = x_mat.ncol();

  const double *xptr   = x_mat.begin();
  const double *v0ptr  = v0_vec.begin();
  const double *v1ptr  = v1_vec.begin();
  const double *etaptr = eta_vec.begin();
  const double *pptr   = params.begin();

  const double beta0 = pptr[0];
  const double betaV = pptr[1];

  double sumQ = 0.0;

  for (int m = 0; m < n; ++m) {
    const int i = (idx == nullptr) ? m : (*idx)[m];
    double xb = 0.0;

    for (int j = 0; j < d; ++j) {
      xb += xptr[i + j * N] * pptr[j + 2];
    }

    const double lp0 = beta0 + betaV * v0ptr[i] + xb;
    const double lp1 = beta0 + betaV * v1ptr[i] + xb;

    const double f0 = mean_function(lp0, family);
    const double f1 = mean_function(lp1, family);

    const double lower = std::min(f0, f1);
    const double upper = std::max(f0, f1);

    const double eta = etaptr[i];

    if (eta < lower) {
      const double diff = lower - eta;
      sumQ += diff * diff;
    } else if (eta > upper) {
      const double diff = eta - upper;
      sumQ += diff * diff;
    }
  }

  return sumQ / static_cast<double>(n);
}

// [[Rcpp::export]]
double Q_obj_cpp(const NumericMatrix &x_mat,
                 const NumericVector &v0_vec,
                 const NumericVector &v1_vec,
                 const NumericVector &params,
                 const NumericVector &eta_vec,
                 const int family = 0) {
  return compute_objective(
    x_mat, v0_vec, v1_vec, params, eta_vec, nullptr,
    static_cast<ModelFamily>(family)
  );
}

// [[Rcpp::export]]
double Q_obj_subsample_cpp(const NumericMatrix &x_mat,
                           const NumericVector &v0_vec,
                           const NumericVector &v1_vec,
                           const NumericVector &params,
                           const NumericVector &eta_vec,
                           const IntegerVector &idx_vec,
                           const int family = 0) {
  return compute_objective(
    x_mat, v0_vec, v1_vec, params, eta_vec, &idx_vec,
    static_cast<ModelFamily>(family)
  );
}

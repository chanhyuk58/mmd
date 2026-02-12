#include <Rcpp.h>
using namespace Rcpp;

// ------------------ Q objective ------------------
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

  double gamma0 = pptr[0];
  double gammaV = pptr[1];
  double sumQ = 0.0;

  for (int i = 0; i < n; ++i) {
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

  return sumQ / static_cast<double>(n);
}

// ------------------ gradient (central FD) ------------------
// [[Rcpp::export]]
NumericVector Q_grad_cpp(const NumericMatrix & x_mat,
                         const NumericVector & v0_vec,
                         const NumericVector & v1_vec,
                         const NumericVector & params,
                         const NumericVector & eta_vec,
                         double rel_step = 1e-6) {

  int p = params.size();
  NumericVector grad(p);
  NumericVector par_plus = clone(params);
  NumericVector par_minus = clone(params);

  for (int k = 0; k < p; ++k) {
    double h = rel_step * std::max(1.0, std::abs((double)params[k]));
    if (h == 0.0) h = rel_step;

    par_plus[k] = params[k] + h;
    par_minus[k] = params[k] - h;

    double q_plus = Q_obj_cpp(x_mat, v0_vec, v1_vec, par_plus, eta_vec);
    double q_minus = Q_obj_cpp(x_mat, v0_vec, v1_vec, par_minus, eta_vec);

    grad[k] = (q_plus - q_minus) / (2.0 * h);

    // restore
    par_plus[k] = params[k];
    par_minus[k] = params[k];
  }

  return grad;
}

// ------------------ simple BFGS optimizer with Armijo line search ------------------
// Minimizes f(par) using gradient g(par) provided via Q_grad_cpp and Q_obj_cpp
// [[Rcpp::export]]
List bfgs_opt_cpp(const NumericMatrix & x_mat,
                  const NumericVector & v0_vec,
                  const NumericVector & v1_vec,
                  const NumericVector & eta_vec,
                  NumericVector par0,
                  int maxiter = 1000,
                  double gtol = 1e-6,
                  double ftol = 1e-12,
                  double wolfe_c1 = 1e-4,
                  double alpha_init = 1.0) {

  int p = par0.size();
  NumericVector x = clone(par0);
  NumericVector grad = Q_grad_cpp(x_mat, v0_vec, v1_vec, x, eta_vec);
  double fval = Q_obj_cpp(x_mat, v0_vec, v1_vec, x, eta_vec);

  // initialize inverse Hessian approximation as identity
  NumericMatrix H(p, p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      H(i, j) = (i == j) ? 1.0 : 0.0;
    }
  }

  int iter = 0;
  double grad_norm = 0.0;
  for (int k = 0; k < p; ++k) grad_norm += grad[k] * grad[k];
  grad_norm = sqrt(grad_norm);

  while (iter < maxiter && grad_norm > gtol) {
    // compute search direction: p = -H * grad
    NumericVector pdir(p);
    for (int i = 0; i < p; ++i) {
      double s = 0.0;
      for (int j = 0; j < p; ++j) s += H(i, j) * grad[j];
      pdir[i] = -s;
    }

    // line search (backtracking Armijo)
    double alpha = alpha_init;
    double f0 = fval;
    double gtp = 0.0;
    for (int j = 0; j < p; ++j) gtp += grad[j] * pdir[j];

    // if gtp >= 0, take steepest descent direction instead
    if (gtp >= 0) {
      for (int j = 0; j < p; ++j) pdir[j] = -grad[j];
      gtp = 0.0;
      for (int j = 0; j < p; ++j) gtp += grad[j] * pdir[j];
    }

    bool armijo_ok = false;
    int ls_iter = 0;
    while (ls_iter < 60) {
      NumericVector xnew = clone(x);
      for (int j = 0; j < p; ++j) xnew[j] = x[j] + alpha * pdir[j];

      double fnew = Q_obj_cpp(x_mat, v0_vec, v1_vec, xnew, eta_vec);

      if (fnew <= f0 + wolfe_c1 * alpha * gtp) { // Armijo condition
        // accept
        // compute new grad
        NumericVector gradnew = Q_grad_cpp(x_mat, v0_vec, v1_vec, xnew, eta_vec);

        // BFGS update
        NumericVector svec(p), yvec(p);
        for (int j = 0; j < p; ++j) {
          svec[j] = alpha * pdir[j];
          yvec[j] = gradnew[j] - grad[j];
        }

        double sty = 0.0;
        for (int j = 0; j < p; ++j) sty += svec[j] * yvec[j];

        if (sty > 1e-12) {
          // H <- (I - rho*s*y^T) H (I - rho*y*s^T) + rho * s * s^T
          double rho = 1.0 / sty;

          // compute (I - rho*s*y^T)
          NumericMatrix Iminus(p, p);
          for (int a = 0; a < p; ++a) {
            for (int b = 0; b < p; ++b) {
              double tmp = 0.0;
              if (a == b) tmp = 1.0;
              tmp -= rho * svec[a] * yvec[b];
              Iminus(a, b) = tmp;
            }
          }

          // temp = Iminus * H
          NumericMatrix temp(p, p);
          for (int a = 0; a < p; ++a) {
            for (int b = 0; b < p; ++b) {
              double sum = 0.0;
              for (int c = 0; c < p; ++c) sum += Iminus(a, c) * H(c, b);
              temp(a, b) = sum;
            }
          }

          // H_new = temp * transpose(Iminus)
          NumericMatrix Hnew(p, p);
          for (int a = 0; a < p; ++a) {
            for (int b = 0; b < p; ++b) {
              double sum = 0.0;
              for (int c = 0; c < p; ++c) sum += temp(a, c) * Iminus(b, c); // transposed index
              Hnew(a, b) = sum;
            }
          }

          // add rho * s s^T
          for (int a = 0; a < p; ++a) {
            for (int b = 0; b < p; ++b) {
              Hnew(a, b) += rho * svec[a] * svec[b];
            }
          }

          H = Hnew;
        } else {
          // skip update if sty is too small; reset H to identity to maintain positive-definiteness
          for (int a = 0; a < p; ++a) {
            for (int b = 0; b < p; ++b) {
              H(a, b) = (a == b) ? 1.0 : 0.0;
            }
          }
        }

        // accept step
        x = clone(xnew);
        grad = clone(gradnew);
        fval = fnew;
        armijo_ok = true;
        break;
      } else {
        alpha *= 0.5; // reduce step
        ls_iter++;
      }
    } // end line search

    if (!armijo_ok) {
      // If line search failed, stop
      break;
    }

    // update gradient norm and iteration
    grad_norm = 0.0;
    for (int j = 0; j < p; ++j) grad_norm += grad[j] * grad[j];
    grad_norm = sqrt(grad_norm);

    if (std::abs(fval) < ftol) break;

    iter++;
  } // end main loop

  return List::create(Named("par") = x,
                      Named("value") = fval,
                      Named("grad") = grad,
                      Named("iter") = iter,
                      Named("conv") = (grad_norm <= gtol));
}

// ------------------ coordinatewise bounds via bisection ------------------
// For coordinate j, we find t in [lo, hi] such that Q(par with par[j]=t) = threshold.
// If the whole interval is feasible (Q<=threshold), we return the endpoint lo/hi.
// [[Rcpp::export]]
List coord_bounds_bisect_cpp(const NumericMatrix & x_mat,
                             const NumericVector & v0_vec,
                             const NumericVector & v1_vec,
                             const NumericVector & eta_vec,
                             const NumericVector & par_opt,
                             double threshold,
                             double box_radius = 1.0,
                             int max_bisect_iter = 60,
                             double tol = 1e-8) {

  int p = par_opt.size();
  NumericVector lower(p), upper(p);

  for (int j = 0; j < p; ++j) {
    // define interval for lower bound: [par_opt[j] - box, par_opt[j]]
    double left = par_opt[j] - box_radius;
    double right = par_opt[j];
    NumericVector parL = clone(par_opt);
    NumericVector parR = clone(par_opt);

    // Evaluate at endpoints
    parL[j] = left;
    parR[j] = right;
    double qL = Q_obj_cpp(x_mat, v0_vec, v1_vec, parL, eta_vec) - threshold;
    double qR = Q_obj_cpp(x_mat, v0_vec, v1_vec, parR, eta_vec) - threshold;

    if (qL <= 0.0) {
      // entire interval feasible on left side -> lower bound = left
      lower[j] = left;
    } else if (qR > 0.0) {
      // weird: even at par_opt (right) Q > threshold (should not happen) -> NA
      lower[j] = NA_REAL;
    } else {
      // qL > 0, qR <= 0, so there is a root in (left, right]
      double a = left;
      double b = right;
      double fa = qL;
      double fb = qR;
      double mid = a;
      for (int it = 0; it < max_bisect_iter; ++it) {
        mid = 0.5 * (a + b);
        NumericVector parMid = clone(par_opt);
        parMid[j] = mid;
        double fm = Q_obj_cpp(x_mat, v0_vec, v1_vec, parMid, eta_vec) - threshold;
        if (std::abs(fm) < tol) break;
        if (fm > 0) {
          a = mid;
          fa = fm;
        } else {
          b = mid;
          fb = fm;
        }
      }
      lower[j] = mid;
    }

    // Now upper bound: interval [par_opt[j], par_opt[j] + box_radius]
    left = par_opt[j];
    right = par_opt[j] + box_radius;
    parL = clone(par_opt); parR = clone(par_opt);
    parL[j] = left; parR[j] = right;
    qL = Q_obj_cpp(x_mat, v0_vec, v1_vec, parL, eta_vec) - threshold;
    qR = Q_obj_cpp(x_mat, v0_vec, v1_vec, parR, eta_vec) - threshold;

    if (qR <= 0.0) {
      // entire interval feasible on right side -> upper bound = right
      upper[j] = right;
    } else if (qL > 0.0) {
      // even at par_opt Q > threshold -> NA
      upper[j] = NA_REAL;
    } else {
      // qL <= 0, qR > 0 => root in (left, right]
      double a = left;
      double b = right;
      double mid = a;
      for (int it = 0; it < max_bisect_iter; ++it) {
        mid = 0.5 * (a + b);
        NumericVector parMid = clone(par_opt);
        parMid[j] = mid;
        double fm = Q_obj_cpp(x_mat, v0_vec, v1_vec, parMid, eta_vec) - threshold;
        if (std::abs(fm) < tol) break;
        if (fm > 0) {
          b = mid;
        } else {
          a = mid;
        }
      }
      upper[j] = mid;
    }
  } // end for j

  return List::create(Named("lower") = lower, Named("upper") = upper);
}

##############################################################################
# hw2_problem2_run.R
# Homework 2 - Problem 2: 2-dim Numerical Integration
# Run on Hellbender via: Rscript hw2_problem2_run.R
##############################################################################

library(Rcpp)

cat("=== HW2 Problem 2: 2-dim Numerical Integration ===\n")
cat("Running on:", Sys.info()["nodename"], "\n")
cat("Date:", as.character(Sys.time()), "\n\n")

# ============================================================
# All C++ in one block:
#   - Incomplete gamma via adaptive GL quadrature
#     (substitution u = t^s for s < 1 singularity)
#   - Upper gamma for negative s via recurrence
#   - Survival function S*(x) and CDF F*(x)
#   - Convolution integral F_Y(y) with survival trick
# ============================================================

sourceCpp(code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// ====== Gauss-Legendre 15-point quadrature ======
static const int GL_N = 15;
static const double gl_nodes[15] = {
  -0.987992518020485, -0.937273392400706, -0.848206583410427,
  -0.724417731360170, -0.570972172608539, -0.394151347077563,
  -0.201194093997435,  0.000000000000000,  0.201194093997435,
   0.394151347077563,  0.570972172608539,  0.724417731360170,
   0.848206583410427,  0.937273392400706,  0.987992518020485
};
static const double gl_weights[15] = {
   0.030753241996117,  0.070366047488108,  0.107159220467172,
   0.139570677926154,  0.166269205816994,  0.186161000015562,
   0.198431485327111,  0.202578241925561,  0.198431485327111,
   0.186161000015562,  0.166269205816994,  0.139570677926154,
   0.107159220467172,  0.070366047488108,  0.030753241996117
};

// ---- Original integrand: t^{s-1} * e^{-t} in log-space ----
inline double integrand_orig(double t, double s) {
  if (t <= 0.0) return 0.0;
  return exp((s - 1.0) * log(t) - t);
}

// ---- Transformed integrand for s < 1 (substitution u = t^s) ----
// After substitution: (1/s) * exp(-u^{1/s})
inline double integrand_transformed(double u, double s) {
  if (u <= 0.0) return 0.0;
  double inv_s = 1.0 / s;
  return inv_s * exp(-pow(u, inv_s));
}

// ---- Function pointer type ----
typedef double (*integrand_fn)(double, double);

// ---- GL quadrature on [a, b] ----
double gl_quad(double a, double b, double s, integrand_fn f) {
  double mid = 0.5 * (a + b);
  double half = 0.5 * (b - a);
  double sum = 0.0;
  for (int i = 0; i < GL_N; i++) {
    double t = mid + half * gl_nodes[i];
    sum += gl_weights[i] * f(t, s);
  }
  return half * sum;
}

// ---- Adaptive quadrature with relative tolerance ----
double adaptive_quad(double a, double b, double s, integrand_fn f,
                     double tol, int depth) {
  double whole = gl_quad(a, b, s, f);
  double mid_pt = 0.5 * (a + b);
  double left_val = gl_quad(a, mid_pt, s, f);
  double right_val = gl_quad(mid_pt, b, s, f);
  double refined = left_val + right_val;
  double err = fabs(refined - whole);
  if (err < tol * (1.0 + fabs(refined)) || depth >= 25) {
    return refined;
  }
  return adaptive_quad(a, mid_pt, s, f, tol, depth + 1) +
         adaptive_quad(mid_pt, b, s, f, tol, depth + 1);
}

// ====== Lower incomplete gamma ======
// For s < 1: split at min(1, x).
//   [0, split]: substitution u = t^s removes singularity
//   [split, x]: direct quadrature (smooth)
// For s >= 1: direct quadrature on [0, x].
// For large x (> max(50, s+10*sqrt(s))): use Gamma(s) - tiny tail.
double lower_gamma_single(double s, double x) {
  if (x <= 0.0) return 0.0;

  // For large x, gamma(s,x) ≈ Gamma(s). Compute via Gamma(s) - tail.
  if (s > 0.0 && x > std::max(50.0, s + 10.0 * sqrt(s))) {
    double tail_len = std::max(50.0, 4.0 * sqrt(s > 1.0 ? s : 1.0));
    double upper_tail = adaptive_quad(x, x + tail_len, s, integrand_orig, 1e-10, 0);
    return R::gammafn(s) - upper_tail;
  }

  if (s < 1.0) {
    double split = std::min(1.0, x);
    double upper_u = pow(split, s);
    double part1 = adaptive_quad(0.0, upper_u, s, integrand_transformed, 1e-10, 0);
    double part2 = 0.0;
    if (x > split) {
      part2 = adaptive_quad(split, x, s, integrand_orig, 1e-10, 0);
    }
    return part1 + part2;
  } else {
    return adaptive_quad(0.0, x, s, integrand_orig, 1e-10, 0);
  }
}

// ====== Upper incomplete gamma — supports negative s ======
// s > 0:  Gamma(s,x) = Gamma(s) - gamma(s,x)
//         For very large x, Gamma(s,x) ≈ 0 (all mass below x).
// s = 0:  direct integration
// s < 0:  recurrence: Gamma(s,x) = (Gamma(s+1,x) - x^s * e^{-x}) / s
double upper_gamma_single(double s, double x) {
  if (x <= 0.0) {
    if (s > 0.0) return R::gammafn(s);
    return R_PosInf;
  }
  // For very large x with s > 0, upper gamma is essentially 0
  if (s > 0.0 && x > std::max(50.0, s + 10.0 * sqrt(s))) {
    // Direct quadrature of the tiny tail
    double tail_len = std::max(50.0, 4.0 * sqrt(s > 1.0 ? s : 1.0));
    return adaptive_quad(x, x + tail_len, s, integrand_orig, 1e-10, 0);
  }
  if (s > 0.0) {
    return R::gammafn(s) - lower_gamma_single(s, x);
  }
  if (fabs(s) < 1e-14) {
    double upper_lim = x + std::max(100.0, 3.0 * x);
    return adaptive_quad(x, upper_lim, 0.0, integrand_orig, 1e-10, 0);
  }
  // s < 0: recurrence Gamma(s,x) = (Gamma(s+1,x) - x^s * e^{-x}) / s
  double upper_sp1 = upper_gamma_single(s + 1.0, x);
  // For very large x, x^s * exp(-x) underflows to 0 when s < 0
  double correction = 0.0;
  if (x < 700.0) {  // exp(-700) ≈ 0, avoid underflow
    correction = pow(x, s) * exp(-x);
  }
  return (upper_sp1 - correction) / s;
}

// [[Rcpp::export]]
NumericVector lower_gamma_cpp(NumericVector s, NumericVector x) {
  int n = s.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) result[i] = lower_gamma_single(s[i], x[i]);
  return result;
}

// [[Rcpp::export]]
NumericVector upper_gamma_cpp(NumericVector s, NumericVector x) {
  int n = s.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) result[i] = upper_gamma_single(s[i], x[i]);
  return result;
}

// ====== Survival function S*(x) = 1 - F*(x) ======
double survival_single(double x, double phi, double tau) {
  if (x <= 0.0) return 1.0;
  double r1pi = 1.0 / sqrt(M_PI);
  double tau2 = tau * tau;
  double arg = pow(tau2 / x, 1.0 / phi);
  double term1 = r1pi * lower_gamma_single(0.5, arg);
  double term2 = (1.0 / x) * r1pi * (tau2 / phi) * upper_gamma_single(0.5 - phi, arg);
  double surv = term1 + term2;
  // Clamp to [0, 1] — small quadrature errors can push slightly outside
  if (surv < 0.0) surv = 0.0;
  if (surv > 1.0) surv = 1.0;
  return surv;
}

double cdf_single(double x, double phi, double tau) {
  if (x <= 0.0) return 0.0;
  return 1.0 - survival_single(x, phi, tau);
}

// [[Rcpp::export]]
NumericVector survival_Fstar(NumericVector x, double phi, double tau) {
  int n = x.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) result[i] = survival_single(x[i], phi, tau);
  return result;
}

// [[Rcpp::export]]
NumericVector cdf_Fstar(NumericVector x, double phi, double tau) {
  int n = x.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) result[i] = cdf_single(x[i], phi, tau);
  return result;
}

// ====== Convolution integral F_Y(y) ======
// F_Y(y) = int_0^inf F*(x) * (1/sigma) * phi((y-x)/sigma) dx
//
// For large y where F*(x)~1: compute survival form instead
// 1 - F_Y(y) = int_0^inf S*(x) * phi((y-x)/sigma)/sigma dx

// [[Rcpp::export]]
NumericVector conv_FY(NumericVector y, double phi, double tau, double sigma,
                      double K_trunc = 8.0) {
  int ny = y.size();
  NumericVector result(ny);
  double inv_sqrt2pi = 1.0 / sqrt(2.0 * M_PI);

  for (int j = 0; j < ny; j++) {
    double yj = y[j];

    double lo = std::max(1e-10, yj - K_trunc * sigma);
    double hi = yj + K_trunc * sigma;

    if (hi <= 1e-10) {
      result[j] = 0.0;
      continue;
    }
    if (lo < 1e-10) lo = 1e-10;

    bool use_survival = (cdf_single(lo, phi, tau) > 0.99);

    double width = hi - lo;
    int n_sub = std::max(1, std::min(1000, (int)ceil(width / sigma)));
    double sub_w = width / (double)n_sub;

    double total = 0.0;
    for (int isub = 0; isub < n_sub; isub++) {
      double a = lo + isub * sub_w;
      double b = a + sub_w;
      double mid = 0.5 * (a + b);
      double half = 0.5 * (b - a);
      double sub_sum = 0.0;

      for (int ig = 0; ig < GL_N; ig++) {
        double xi = mid + half * gl_nodes[ig];
        if (xi <= 0.0) continue;

        double z = (yj - xi) / sigma;
        double kern = inv_sqrt2pi * exp(-0.5 * z * z) / sigma;

        if (use_survival) {
          sub_sum += gl_weights[ig] * survival_single(xi, phi, tau) * kern;
        } else {
          sub_sum += gl_weights[ig] * cdf_single(xi, phi, tau) * kern;
        }
      }
      total += half * sub_sum;
    }

    if (use_survival) {
      result[j] = 1.0 - total;
    } else {
      result[j] = total;
    }

    if (result[j] < 0.0) result[j] = 0.0;
    if (result[j] > 1.0) result[j] = 1.0;
  }
  return result;
}
')

cat("C++ functions compiled successfully.\n\n")

# ============================================================
# Task 1: Verify F*(x) is a valid CDF for (phi, tau) = (0.7, 2)
# ============================================================

cat("=== Task 1: Verify F*(x) is a valid CDF ===\n")
cat("Parameters: phi = 0.7, tau = 2\n\n")

phi_val <- 0.7
tau_val <- 2.0

x_grid <- c(seq(1e-6, 1, length.out = 500),
            seq(1, 100, length.out = 500),
            seq(100, 10000, length.out = 500),
            seq(10000, 1e6, length.out = 200))

Fstar_vals <- cdf_Fstar(x_grid, phi_val, tau_val)

cat("(1) Range check: 0 <= F*(x) <= 1\n")
cat("    min(F*):", min(Fstar_vals), "\n")
cat("    max(F*):", max(Fstar_vals), "\n")
cat("    All in [0,1]?", all(Fstar_vals >= -1e-10 & Fstar_vals <= 1 + 1e-10), "\n\n")

diffs <- diff(Fstar_vals)
cat("(2) Non-decreasing check:\n")
cat("    min(diff):", min(diffs), "\n")
cat("    Any decrease > 1e-10?", any(diffs < -1e-10), "\n")
if (any(diffs < -1e-10)) {
  bad_idx <- which(diffs < -1e-10)
  cat("    Violations at indices:", head(bad_idx, 10), "\n")
  cat("    x values:", head(x_grid[bad_idx], 10), "\n")
  cat("    F* values:", head(Fstar_vals[bad_idx], 10), "\n")
  cat("    F* next:  ", head(Fstar_vals[bad_idx + 1], 10), "\n")
}
cat("\n")

cat("(3) Limit checks:\n")
Fstar_near0 <- cdf_Fstar(c(1e-10, 1e-8, 1e-6), phi_val, tau_val)
cat("    F*(1e-10) =", Fstar_near0[1], " (should -> 0)\n")
cat("    F*(1e-8)  =", Fstar_near0[2], "\n")
cat("    F*(1e-6)  =", Fstar_near0[3], "\n")

Fstar_large <- cdf_Fstar(c(1e4, 1e5, 1e6, 1e8), phi_val, tau_val)
cat("    F*(1e4)   =", Fstar_large[1], "\n")
cat("    F*(1e5)   =", Fstar_large[2], "\n")
cat("    F*(1e6)   =", Fstar_large[3], "\n")
cat("    F*(1e8)   =", Fstar_large[4], " (should -> 1)\n\n")

cat("Negative s check (s = 0.5 - 0.7 = -0.2):\n")
cat("  Gamma(-0.2, 1.0) =", upper_gamma_cpp(-0.2, 1.0), "\n")
cat("  Gamma(-0.2, 5.0) =", upper_gamma_cpp(-0.2, 5.0), "\n\n")

# ============================================================
# Task 1: Plot F*(x)
# ============================================================

pdf("hw2_prob2_task1_CDF_check.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

x1 <- seq(0.001, 100, length.out = 1000)
plot(x1, cdf_Fstar(x1, 0.7, 2), type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(F^"*"*(x)),
     main = expression(paste(phi, " = 0.7, ", tau, " = 2")))
abline(h = c(0, 1), lty = 2, col = "gray")

x2 <- seq(1, 10000, length.out = 1000)
plot(x2, cdf_Fstar(x2, 0.7, 2), type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(F^"*"*(x)),
     main = expression(paste(phi, " = 0.7, ", tau, " = 2, heavy tail")))
abline(h = 1, lty = 2, col = "gray")

dev.off()
cat("Saved: hw2_prob2_task1_CDF_check.pdf\n\n")

# ============================================================
# Task 2: Convolution integral F_Y(y)
# ============================================================

cat("=== Task 2: Convolution CDF F_Y(y) ===\n\n")

y_vals1 <- seq(-100, 100, length.out = 2000)
y_vals2 <- seq(1000, 10000, length.out = 2000)

# Setting 1: (phi, tau, sigma) = (0.3, 2, 3)
cat("--- (phi, tau, sigma) = (0.3, 2, 3) ---\n"); flush.console()
ptm <- proc.time()
FY_03_y1 <- conv_FY(y_vals1, phi = 0.3, tau = 2, sigma = 3, K_trunc = 8)
FY_03_y2 <- conv_FY(y_vals2, phi = 0.3, tau = 2, sigma = 3, K_trunc = 8)
time_03 <- (proc.time() - ptm)[3]
cat("Wall-clock time (4000 points):", time_03, "sec\n")

# Setting 2: (phi, tau, sigma) = (0.7, 2, 3)
cat("--- (phi, tau, sigma) = (0.7, 2, 3) ---\n"); flush.console()
ptm <- proc.time()
FY_07_y1 <- conv_FY(y_vals1, phi = 0.7, tau = 2, sigma = 3, K_trunc = 8)
FY_07_y2 <- conv_FY(y_vals2, phi = 0.7, tau = 2, sigma = 3, K_trunc = 8)
time_07 <- (proc.time() - ptm)[3]
cat("Wall-clock time (4000 points):", time_07, "sec\n\n")

# ============================================================
# Task 2: Reproduce Figure 1 (four panels)
# ============================================================

pdf("hw2_prob2_task2_FY_CDF.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

plot(y_vals1, FY_03_y1, type = "l", col = "blue", lwd = 2,
     xlab = "y", ylab = expression(F[Y](y)),
     main = expression(paste("(a) ", (phi*", "*tau*", "*sigma), " = (0.3, 2, 3), y " %in% " [-100, 100]")))
abline(h = c(0, 1), lty = 2, col = "gray")

plot(y_vals1, FY_07_y1, type = "l", col = "blue", lwd = 2,
     xlab = "y", ylab = expression(F[Y](y)),
     main = expression(paste("(b) ", (phi*", "*tau*", "*sigma), " = (0.7, 2, 3), y " %in% " [-100, 100]")))
abline(h = c(0, 1), lty = 2, col = "gray")

plot(y_vals2, FY_03_y2, type = "l", col = "blue", lwd = 2,
     xlab = "y", ylab = expression(F[Y](y)),
     main = expression(paste("(c) ", (phi*", "*tau*", "*sigma), " = (0.3, 2, 3), y " %in% " [10"^3*", 10"^4*"]")))
abline(h = 1, lty = 2, col = "gray")

plot(y_vals2, FY_07_y2, type = "l", col = "blue", lwd = 2,
     xlab = "y", ylab = expression(F[Y](y)),
     main = expression(paste("(d) ", (phi*", "*tau*", "*sigma), " = (0.7, 2, 3), y " %in% " [10"^3*", 10"^4*"]")))
abline(h = 1, lty = 2, col = "gray")

dev.off()
cat("Saved: hw2_prob2_task2_FY_CDF.pdf\n\n")

# ============================================================
# Monotonicity checks
# ============================================================

cat("=== Monotonicity checks on F_Y ===\n")
cat("phi=0.3, y in [-100,100]:  min(diff) =", min(diff(FY_03_y1)), "\n")
cat("phi=0.3, y in [1e3, 1e4]:  min(diff) =", min(diff(FY_03_y2)), "\n")
cat("phi=0.7, y in [-100,100]:  min(diff) =", min(diff(FY_07_y1)), "\n")
cat("phi=0.7, y in [1e3, 1e4]:  min(diff) =", min(diff(FY_07_y2)), "\n\n")

# ============================================================
# Runtime summary
# ============================================================

cat("=== Runtime Summary ===\n")
cat("(phi=0.3): ", time_03, "sec for 4000 F_Y evaluations\n")
cat("(phi=0.7): ", time_07, "sec for 4000 F_Y evaluations\n")

cat("\nDone!\n")

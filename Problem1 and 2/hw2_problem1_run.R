##############################################################################
# hw2_problem1_run.R
# Homework 2 - Problem 1: Incomplete Gamma Functions
# Run this on Hellbender via: Rscript hw2_problem1_run.R
##############################################################################

library(Rcpp)

cat("=== HW2 Problem 1: Incomplete Gamma Functions ===\n")
cat("Running on:", Sys.info()["nodename"], "\n")
cat("Date:", as.character(Sys.time()), "\n\n")

# ============================================================
# Part (1): Pure R implementation
# ============================================================

gamma_integrand <- function(t, s) {
  t^(s - 1) * exp(-t)
}

lower_gamma_r <- function(s, x) {
  mapply(function(si, xi) {
    if (xi <= 0) return(0)
    tryCatch(
      integrate(gamma_integrand, lower = 0, upper = xi, s = si,
                rel.tol = 1e-10, abs.tol = 1e-14,
                subdivisions = 500L)$value,
      error = function(e) NA_real_
    )
  }, s, x)
}

upper_gamma_r <- function(s, x) {
  mapply(function(si, xi) {
    if (xi <= 0) return(gamma(si))
    tryCatch(
      integrate(gamma_integrand, lower = xi, upper = Inf, s = si,
                rel.tol = 1e-10, abs.tol = 1e-14,
                subdivisions = 500L)$value,
      error = function(e) {
        gamma(si) - lower_gamma_r(si, xi)
      }
    )
  }, s, x)
}

# ============================================================
# Part (2): C++ integrand, R integration
# ============================================================

cppFunction('
NumericVector gamma_integrand_cpp(NumericVector t, double s) {
  int n = t.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    if (t[i] <= 0.0) {
      result[i] = (s >= 1.0) ? 0.0 : R_PosInf;
    } else {
      result[i] = exp((s - 1.0) * log(t[i]) - t[i]);
    }
  }
  return result;
}
')

lower_gamma_rcpp <- function(s, x) {
  mapply(function(si, xi) {
    if (xi <= 0) return(0)
    f <- function(t) gamma_integrand_cpp(t, si)
    tryCatch(
      integrate(f, lower = 0, upper = xi,
                rel.tol = 1e-10, abs.tol = 1e-14,
                subdivisions = 500L)$value,
      error = function(e) NA_real_
    )
  }, s, x)
}

upper_gamma_rcpp <- function(s, x) {
  mapply(function(si, xi) {
    if (xi <= 0) return(gamma(si))
    f <- function(t) gamma_integrand_cpp(t, si)
    tryCatch(
      integrate(f, lower = xi, upper = Inf,
                rel.tol = 1e-10, abs.tol = 1e-14,
                subdivisions = 500L)$value,
      error = function(e) {
        gamma(si) - lower_gamma_rcpp(si, xi)
      }
    )
  }, s, x)
}

# ============================================================
# Part (3): BONUS — Both integrand and numerical integration in C++
#
# Strategy:
#   - Integrand: t^{s-1} * e^{-t} computed in log-space for stability.
#   - Numerical integration: adaptive Gauss-Legendre 15-point quadrature.
#   - Singularity handling (s < 1): the integrand blows up at t = 0.
#     We apply the substitution u = t^s, so t = u^{1/s}, dt = (1/s)*u^{1/s - 1} du.
#     The integral becomes: integral_0^{x^s} (1/s) * e^{-u^{1/s}} du,
#     which is bounded and smooth on [0, x^s].
#   - For the upper incomplete gamma, we integrate from x to a finite cutoff
#     (tail of t^{s-1} e^{-t} decays exponentially).
# ============================================================

sourceCpp(code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Gauss-Legendre 15-point nodes and weights on [-1, 1]
// Roots of the 15th Legendre polynomial; exact for poly degree <= 29.
// Standard mathematical constants (Abramowitz & Stegun Table 25.4).
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
inline double integrand(double t, double s) {
  if (t <= 0.0) return 0.0;
  return exp((s - 1.0) * log(t) - t);
}

// ---- Transformed integrand for s < 1 (substitution u = t^s) ----
// After substitution: (1/s) * exp(-u^{1/s})
// This removes the singularity at t = 0.
inline double integrand_transformed(double u, double s) {
  if (u <= 0.0) return 0.0;
  double inv_s = 1.0 / s;
  return inv_s * exp(-pow(u, inv_s));
}

// ---- Generic GL quadrature on [a, b] with a function pointer ----
typedef double (*integrand_fn)(double, double);

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

// ---- Lower incomplete gamma: integral_0^x t^{s-1} e^{-t} dt ----
// For s < 1, the integrand has a singularity at t = 0.
// Strategy: split into [0, min(1,x)] and [min(1,x), x].
//   - On [0, min(1,x)]: apply substitution u = t^s to remove singularity.
//     Integral becomes: integral_0^{min(1,x)^s} (1/s) exp(-u^{1/s}) du
//   - On [min(1,x), x]: integrand is smooth, use direct quadrature.
// For s >= 1: no singularity, integrate directly on [0, x].
// For large x (> 50): gamma(s,x) ≈ Gamma(s) since tail is negligible.
//   Compute as Gamma(s) - integral_x^{x+tail} t^{s-1} e^{-t} dt.
double lower_gamma_single(double s, double x) {
  if (x <= 0.0) return 0.0;

  // For large x, the integral is essentially Gamma(s).
  // Compute via: Gamma(s) - upper tail (which is tiny).
  if (s > 0.0 && x > std::max(50.0, s + 10.0 * sqrt(s))) {
    // The mode of t^{s-1}e^{-t} is at t = s-1 (for s>1) or near 0 (for s<1).
    // For x >> mode, the tail integral from x to inf is tiny.
    double tail_len = std::max(50.0, 4.0 * sqrt(s > 1.0 ? s : 1.0));
    double upper_tail = adaptive_quad(x, x + tail_len, s, integrand, 1e-10, 0);
    return R::gammafn(s) - upper_tail;
  }

  if (s < 1.0) {
    double split = std::min(1.0, x);  // split point
    // Part 1: [0, split] via substitution u = t^s
    double upper_u = pow(split, s);
    double part1 = adaptive_quad(0.0, upper_u, s, integrand_transformed, 1e-10, 0);
    // Part 2: [split, x] via direct quadrature (no singularity here)
    double part2 = 0.0;
    if (x > split) {
      part2 = adaptive_quad(split, x, s, integrand, 1e-10, 0);
    }
    return part1 + part2;
  } else {
    // s >= 1: no singularity, integrate directly
    return adaptive_quad(0.0, x, s, integrand, 1e-10, 0);
  }
}

// ---- Upper incomplete gamma: integral_x^inf t^{s-1} e^{-t} dt ----
// Use identity: Gamma(s, x) = Gamma(s) - gamma(s, x)
double upper_gamma_single(double s, double x) {
  if (x <= 0.0) return R::gammafn(s);
  return R::gammafn(s) - lower_gamma_single(s, x);
}

// [[Rcpp::export]]
NumericVector lower_gamma_cpp(NumericVector s, NumericVector x) {
  int n = s.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = lower_gamma_single(s[i], x[i]);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector upper_gamma_cpp(NumericVector s, NumericVector x) {
  int n = s.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = upper_gamma_single(s[i], x[i]);
  }
  return result;
}
')

# ============================================================
# Sanity Checks
# ============================================================

cat("=== Sanity Checks ===\n")
cat("γ(1,1):  R =", lower_gamma_r(1,1), " C++ =", lower_gamma_cpp(1,1),
    " Expected =", 1-exp(-1), "\n")
cat("Γ(1,1):  R =", upper_gamma_r(1,1), " C++ =", upper_gamma_cpp(1,1),
    " Expected =", exp(-1), "\n")
cat("γ(0.5,2): R =", lower_gamma_r(0.5,2), " C++ =", lower_gamma_cpp(0.5,2),
    " Expected =", pgamma(2, 0.5, 1, lower.tail=TRUE)*gamma(0.5), "\n")

# Edge cases that caused the timeout before
cat("\n--- Edge case checks (these caused timeout before) ---\n")
cat("γ(0.0001, 0.0001): C++ =", lower_gamma_cpp(0.0001, 0.0001),
    " Baseline =", pgamma(0.0001, 0.0001, 1, lower.tail=TRUE)*gamma(0.0001), "\n")
cat("γ(0.0001, 30):     C++ =", lower_gamma_cpp(0.0001, 30),
    " Baseline =", pgamma(30, 0.0001, 1, lower.tail=TRUE)*gamma(0.0001), "\n")
cat("γ(10, 30):         C++ =", lower_gamma_cpp(10, 30),
    " Baseline =", pgamma(30, 10, 1, lower.tail=TRUE)*gamma(10), "\n\n")

# ============================================================
# Task 1: Timing over 10,000 (s, x) pairs
# ============================================================

s_vals <- seq(0.0001, 10, length.out = 100)
x_vals <- seq(0.0001, 30, length.out = 100)
param_set <- expand.grid(s = s_vals, x = x_vals)
s_vec <- param_set$s
x_vec <- param_set$x

cat("=== Task 1: Timing (10,000 pairs) ===\n\n")

# --- Method 1: Pure R ---
cat("Running Method 1 (Pure R)...\n"); flush.console()
ptm <- proc.time()
lower_r_results <- lower_gamma_r(s_vec, x_vec)
time_lower_r <- (proc.time() - ptm)[3]

ptm <- proc.time()
upper_r_results <- upper_gamma_r(s_vec, x_vec)
time_upper_r <- (proc.time() - ptm)[3]

cat("Method 1 (Pure R):\n")
cat("  Lower γ(s,x):", time_lower_r, "sec\n")
cat("  Upper Γ(s,x):", time_upper_r, "sec\n\n")

# --- Method 2: Rcpp integrand + R integrate ---
cat("Running Method 2 (Rcpp integrand + R integrate)...\n"); flush.console()
ptm <- proc.time()
lower_rcpp_results <- lower_gamma_rcpp(s_vec, x_vec)
time_lower_rcpp <- (proc.time() - ptm)[3]

ptm <- proc.time()
upper_rcpp_results <- upper_gamma_rcpp(s_vec, x_vec)
time_upper_rcpp <- (proc.time() - ptm)[3]

cat("Method 2 (Rcpp integrand + R integrate):\n")
cat("  Lower γ(s,x):", time_lower_rcpp, "sec\n")
cat("  Upper Γ(s,x):", time_upper_rcpp, "sec\n\n")

# --- Method 3: Full C++ (Bonus) ---
cat("Running Method 3 (Full C++)...\n"); flush.console()
ptm <- proc.time()
lower_cpp_results <- lower_gamma_cpp(s_vec, x_vec)
time_lower_cpp <- (proc.time() - ptm)[3]

ptm <- proc.time()
upper_cpp_results <- upper_gamma_cpp(s_vec, x_vec)
time_upper_cpp <- (proc.time() - ptm)[3]

cat("Method 3 (Full C++ Bonus):\n")
cat("  Lower γ(s,x):", time_lower_cpp, "sec\n")
cat("  Upper Γ(s,x):", time_upper_cpp, "sec\n\n")

# ============================================================
# Task 2: Baseline comparison (MSE)
# ============================================================

cat("=== Task 2: Baseline Comparison ===\n\n")

ptm <- proc.time()
baseline_lower <- pgamma(x_vec, shape = s_vec, scale = 1, lower.tail = TRUE) * gamma(s_vec)
time_bl_lower <- (proc.time() - ptm)[3]

ptm <- proc.time()
baseline_upper <- pgamma(x_vec, shape = s_vec, scale = 1, lower.tail = FALSE) * gamma(s_vec)
time_bl_upper <- (proc.time() - ptm)[3]

cat("Baseline (pgamma * gamma):\n")
cat("  Lower:", time_bl_lower, "sec\n")
cat("  Upper:", time_bl_upper, "sec\n\n")

valid_idx <- is.finite(baseline_lower) & is.finite(baseline_upper)
mse_fn <- function(a, b) mean((a[valid_idx] - b[valid_idx])^2, na.rm = TRUE)

cat("--- MSE vs Baseline ---\n")
cat("Method 1 (Pure R):\n")
cat("  Lower MSE:", mse_fn(lower_r_results, baseline_lower), "\n")
cat("  Upper MSE:", mse_fn(upper_r_results, baseline_upper), "\n")
cat("Method 2 (Rcpp integrand):\n")
cat("  Lower MSE:", mse_fn(lower_rcpp_results, baseline_lower), "\n")
cat("  Upper MSE:", mse_fn(upper_rcpp_results, baseline_upper), "\n")
cat("Method 3 (Full C++):\n")
cat("  Lower MSE:", mse_fn(lower_cpp_results, baseline_lower), "\n")
cat("  Upper MSE:", mse_fn(upper_cpp_results, baseline_upper), "\n\n")

# ============================================================
# Summary Table
# ============================================================

cat("=== Summary Table ===\n")
timing_df <- data.frame(
  Method = c("Pure R", "Rcpp integrand + R integrate",
             "Full C++ (Bonus)", "R Baseline (pgamma)"),
  Lower_sec = c(time_lower_r, time_lower_rcpp, time_lower_cpp, time_bl_lower),
  Upper_sec = c(time_upper_r, time_upper_rcpp, time_upper_cpp, time_bl_upper)
)
print(timing_df, row.names = FALSE)

cat("\nDone!\n")

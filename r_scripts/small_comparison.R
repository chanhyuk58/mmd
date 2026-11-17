# Small comparison based on MMD/draft/note20261016.tex.
# Run from the project root: Rscript --vanilla MMD/r_scripts/small_comparison.R
# Optional first argument: number of repetitions (default 30).
library(splines)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 30L
N <- 100L; TT <- 20L; truth <- 0.08
set.seed(20261023)
out <- file.path("MMD", "output", "small_comparison")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Each country's initial age is uniform on 0:A; 20% enter with known age 0.
# Duration is measured BEFORE drawing the current onset (including event rows).
generate <- function(shape, assignment) {
  A <- sample(10:30, N, replace = TRUE)
  A[sample(N, N / 5)] <- 0
  age <- vapply(A, function(a) sample.int(a + 1, 1) - 1L, integer(1))
  lo <- integer(N); unknown <- A > 0; spell <- integer(N)
  ans <- vector("list", TT)
  for (t in 0:(TT - 1)) {
    z <- runif(N, -1, 1)
    exposure_age <- if (assignment == "observed") A else age
    x <- rbinom(N, 1, plogis(-1 + 3 * exposure_age / 30 + 1.2 * z))
    f <- if (shape == "linear") 0.20 - 0.003 * age else 0.06 + 0.16 * exp(-age / 8)
    p <- f + truth * x + 0.03 * z^2
    stopifnot(all(p > 0 & p < 1))
    y <- rbinom(N, 1, p)
    ans[[t + 1]] <- data.frame(id = 1:N, year = t, y, x, z, age,
      lo, hi = lo + unknown * A, A, unknown, spell = paste(1:N, spell, sep = "_"))
    age <- ifelse(y == 1, 1, age + 1)
    lo <- ifelse(y == 1, 1, lo + 1)
    unknown <- unknown & y == 0
    spell <- spell + y
  }
  do.call(rbind, ans)
}

# Fixed spline basis shared by the oracle, naive, drop, and latent-age fits.
basis <- function(p) bs(p, knots = c(5, 15), Boundary.knots = c(0, 50),
                        degree = 3, intercept = TRUE)
ols <- function(d, p) unname(lm.fit(cbind(x = d$x, z2 = d$z^2, basis(p)), d$y)$coef[1])

# Two folds, split by COUNTRY. Least-squares spline series learn both nuisances.
residuals_cf <- function(response, controls, fold) {
  residual <- matrix(NA_real_, nrow(response), ncol(response))
  for (k in 1:2) {
    test <- fold == k
    fit <- lm.fit(controls[!test, , drop = FALSE], response[!test, , drop = FALSE])
    residual[test, ] <- response[test, , drop = FALSE] - controls[test, , drop = FALSE] %*% fit$coef
  }
  residual
}
partial_out <- function(r) sum(r[, 1] * r[, 2]) / sum(r[, 2]^2)

# Uniform latent-age likelihood: ONE age for the entire initial spell.
# Summing the likelihood over ages differs from averaging fitted coefficients.
latent_spline <- function(d) {
  rows <- groups <- shifts <- owners <- integer()
  offset <- 0L
  for (i in 1:N) {
    idx <- which(d$id == i)
    a <- 0:d$A[idx[1]]
    rows <- c(rows, rep(idx, times = length(a)))
    shifts <- c(shifts, rep(a, each = length(idx)))
    groups <- c(groups, rep(offset + seq_along(a), each = length(idx)))
    owners <- c(owners, rep(i, length(a)))
    offset <- offset + length(a)
  }
  P <- d$lo[rows] + d$unknown[rows] * shifts
  X <- cbind(d$x[rows], d$z[rows]^2, basis(P))
  y <- d$y[rows]
  evaluate <- function(b) {
    p <- as.vector(X %*% b)
    ll <- as.vector(rowsum(y * log(p) + (1 - y) * log1p(-p), groups, reorder = FALSE))
    peak <- ave(ll, owners, FUN = max)
    w <- exp(ll - peak)
    total <- ave(w, owners, FUN = sum)
    post <- w / total
    first <- !duplicated(owners)
    value <- -sum(peak[first] + log(total[first]) - log(tabulate(owners)))
    gradient <- -as.vector(crossprod(X, post[groups] * (y - p) / (p * (1 - p))))
    list(value = value, gradient = gradient)
  }
  # Nonnegative B-splines sum to one; these bounds keep EVERY p in (0,1).
  fit <- optim(c(0.08, 0.03, rep(0.15, 6)),
    function(b) evaluate(b)$value, function(b) evaluate(b)$gradient,
    method = "L-BFGS-B", lower = c(-0.02, -0.01, rep(0.031, 6)),
    upper = c(0.25, 0.10, rep(0.60, 6)), control = list(maxit = 300, factr = 1e7))
  c(estimate = fit$par[1], convergence = fit$convergence)
}

# Linear MMD: profile the interval-violation loss over beta_x.
# The returned range is a criterion set, NOT a 95% confidence interval.
mmd <- function(d, eta) {
  grid <- seq(-0.20, 0.30, by = 0.01)
  objective <- function(b, bx) {
    end0 <- b[1] + bx * d$x + b[2] * d$lo / 30 + b[3] * d$z^2
    end1 <- b[1] + bx * d$x + b[2] * d$hi / 30 + b[3] * d$z^2
    mean(pmax(pmin(end0, end1) - eta, 0)^2 +
         pmax(eta - pmax(end0, end1), 0)^2)
  }
  gradient <- function(b, bx) {
    L <- cbind(1, d$hi / 30, d$z^2); U <- cbind(1, d$lo / 30, d$z^2)
    low <- pmax(as.vector(L %*% b) + bx * d$x - eta, 0)
    high <- pmax(eta - as.vector(U %*% b) - bx * d$x, 0)
    2 * colMeans(L * low - U * high)
  }
  fits <- lapply(grid, function(bx) optim(c(0.20, -0.09, 0.03), objective, gradient, bx = bx,
    method = "L-BFGS-B", lower = c(-1, -1, -1), upper = c(1, 0, 1)))
  Q <- vapply(fits, function(fit) fit$value, numeric(1))
  keep <- Q <= min(Q) + log(nrow(d)) / nrow(d)
  c(lower = min(grid[keep]), upper = max(grid[keep]),
    grid_edge = keep[1] | keep[length(grid)],
    failures = sum(vapply(fits, function(fit) fit$convergence != 0, logical(1))))
}

run_one <- function(shape, assignment, iteration) {
  d <- generate(shape, assignment)
  stopifnot(all(d$lo <= d$age & d$age <= d$hi),
            all(d$age[!d$unknown] == d$lo[!d$unknown]))
  fold <- sample(rep(1:2, length.out = N))[d$id]
  controls <- model.matrix(~ ns(A, 3) + ns(lo, 3) * unknown + ns(z, 3) + year, d)
  r <- residuals_cf(cbind(d$y, d$x), controls, fold)
  # Within-spell transformation stays inside each country's assigned fold.
  W <- cbind(y = d$y, x = d$x,
    model.matrix(~ factor(year) + z + I(z^2) + I(z^3), d)[, -1])
  W <- apply(W, 2, function(v) v - ave(v, d$spell))
  fe <- lm.fit(W[, -1], W[, 1])$coef[1]
  r_fe <- residuals_cf(W[, 1:2], W[, -(1:2)], fold)
  # MMD learns E[Y | observed interval, X, Z] without seeing the true ages.
  eta <- d$y - residuals_cf(matrix(d$y), cbind(controls, x = d$x), fold)[, 1]
  bounds <- mmd(d, pmin(1, pmax(0, eta)))
  latent <- latent_spline(d)
  known <- !d$unknown
  estimates <- c(Oracle = ols(d, d$age), Naive = ols(d, d$lo),
    Drop = ols(d[known, ], d$lo[known]), FE = unname(fe),
    DML = partial_out(r), DML_FE = partial_out(r_fe), Spline_uniform = latent[1])
  names(estimates)[7] <- "Spline_uniform"
  list(points = data.frame(shape, assignment, iteration, method = names(estimates),
                          estimate = unname(estimates)),
    sets = data.frame(shape, assignment, iteration, lower = bounds[1], upper = bounds[2],
      contains_truth = bounds[1] <= truth & truth <= bounds[2],
      grid_edge = bounds[3], optimizer_failures = bounds[4]),
    diagnostics = data.frame(shape, assignment, iteration, onset_rate = mean(d$y),
      unknown_share = mean(d$unknown), spline_convergence = latent[2]))
}

settings <- expand.grid(shape = c("linear", "nonlinear"),
                        assignment = c("observed", "history"), stringsAsFactors = FALSE)
runs <- list()
for (s in seq_len(nrow(settings))) {
  for (b in seq_len(reps)) {
    runs[[length(runs) + 1]] <- run_one(settings$shape[s], settings$assignment[s], b)
    if (b %% 5 == 0 || b == reps)
      cat(settings$shape[s], settings$assignment[s], b, "/", reps, "\n")
  }
}
points <- do.call(rbind, lapply(runs, `[[`, "points"))
sets <- do.call(rbind, lapply(runs, `[[`, "sets"))
diagnostics <- do.call(rbind, lapply(runs, `[[`, "diagnostics"))
summary_points <- do.call(rbind, lapply(split(points, interaction(points$shape, points$assignment, points$method)),
  function(d) data.frame(d[1, c("shape", "assignment", "method")], repetitions = nrow(d),
    mean = mean(d$estimate), bias = mean(d$estimate - truth),
    rmse = sqrt(mean((d$estimate - truth)^2)), mcse_mean = sd(d$estimate) / sqrt(nrow(d)))))
summary_sets <- aggregate(cbind(lower, upper, contains_truth, grid_edge, optimizer_failures) ~
                            shape + assignment, sets, mean)
summary_sets$width <- summary_sets$upper - summary_sets$lower
write.csv(points, file.path(out, "estimates.csv"), row.names = FALSE)
write.csv(summary_points, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(sets, file.path(out, "mmd_sets.csv"), row.names = FALSE)
write.csv(summary_sets, file.path(out, "mmd_summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(summary_points, row.names = FALSE, digits = 3)
print(summary_sets, row.names = FALSE, digits = 3)

png(file.path(out, "comparison.png"), width = 1900, height = 1400, res = 160)
par(mfrow = c(2, 2), mar = c(4.2, 8.8, 3, 1), oma = c(0, 0, 2, 0))
methods <- c("Oracle", "Naive", "Drop", "FE", "DML", "DML_FE", "Spline_uniform")
labels <- c("Oracle", "Zero at entry", "Drop unknown", "Spell + year FE",
            "DML", "DML + FE", "Uniform-age spline")
for (s in seq_len(nrow(settings))) {
  d <- subset(points, shape == settings$shape[s] & assignment == settings$assignment[s])
  boxplot(lapply(methods, function(m) d$estimate[d$method == m]), horizontal = TRUE,
    names = labels, las = 1, col = c("#B8D8BA", rep("#BFD3E6", 6)),
    ylim = range(c(points$estimate, truth)), xlab = "Estimated risk difference",
    main = paste(settings$shape[s], "/", settings$assignment[s]))
  abline(v = truth, col = "#B64242", lty = 2, lwd = 2)
}
mtext("Small Monte Carlo: true effect = 0.08; observed vs. hidden-history assignment",
      outer = TRUE, cex = 1.05)
invisible(dev.off())

# Figures for agent/slides/calibration_opening.tex, including the supplements.
# Run from the project root: Rscript R/calibration_figures.R
# In RStudio, source("R/calibration_figures.R") runs the same sequential script.
# Install once: install.packages(c("mvtnorm", "pracma"))

# 1. Setup ---------------------------------------------------------------

project_dir <- getwd()
figure_dir <- file.path(project_dir, "figures")
result_dir <- file.path(project_dir, "output", "results")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("mvtnorm", quietly = TRUE),
          requireNamespace("pracma", quietly = TRUE), capabilities("cairo"))

# 2. Opening and shared examples -----------------------------------------

message("Opening and shared examples")
cairo_pdf(file.path(figure_dir, "opening_regimes.pdf"), width = 5.1, height = 4,
          family = "sans", pointsize = 16)
par(mar = c(3.2, 3.3, 1.6, 0.7), mgp = c(2.0, 0.6, 0), las = 1,
    bg = "white", fg = "#19324A", col.axis = "#19324A",
    col.lab = "#19324A", col.main = "#19324A",
    family = "sans", xaxs = "i", yaxs = "i", bty = "l")
x <- seq(-4.5, 4.5, length.out = 1201)
left <- dnorm(x, -1, 1)
right <- dnorm(x, 1, 1)
mixture <- (left + right) / 2
plot(x, left, type = "l", lwd = 2.7, col = "#167D8D",
     ylim = c(0, 0.44), xlab = "Outcome y", ylab = "Density",
     axes = FALSE)
axis(1, at = c(-4, -2, 0, 2, 4))
axis(2, at = c(0, 0.2, 0.4), las = 1)
lines(x, right, col = "#C96539", lwd = 2.7)
lines(x, mixture, col = "#626D77", lwd = 2.5, lty = 2)
mtext("Issued forecasts and their mixture", side = 3, line = 0.3,
      cex = 0.94, font = 2)
legend("topright", legend = c("M = -1", "M = +1", "Average"),
       col = c("#167D8D", "#C96539", "#626D77"),
       lty = c(1, 1, 2), lwd = 2.5, bty = "n", cex = 0.80)
box(bty = "l")
dev.off()

cairo_pdf(file.path(figure_dir, "opening_square.pdf"), width = 4.1, height = 3.5,
          family = "sans", pointsize = 17)
t <- .5
stopifnot(length(t) == 1, t > 0, t < 1)
par(mar = c(3.2, 3.4, 0.7, 0.9), mgp = c(2.0, 0.6, 0), las = 1,
    bg = "white", fg = "#19324A", col.axis = "#19324A",
    col.lab = "#19324A", col.main = "#19324A",
    family = "sans", xaxs = "i", yaxs = "i", bty = "l")
par(pty = "s")
plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
     xlab = expression(y[1]), ylab = expression(y[2]), axes = FALSE)
x <- seq(t, 1, length.out = 600)
polygon(c(0, 0, t, x[-1], 1), c(0, 1, 1, t / x[-1], 0),
        col = "#D1E9EA", border = NA)
lines(x, t / x, col = "#167D8D", lwd = 3)
axis(1, at = c(0, 0.5, 1))
axis(2, at = c(0, 0.5, 1), las = 1)
box(bty = "o", col = "#19324A")
area <- t * (1 - log(t))
text(0.34, 0.35, "Shaded area", cex = 0.94)
text(0.34, 0.24, sprintf("%.3f", area), cex = 1.18, font = 2, col = "#167D8D")
text(0.91, 0.88, bquote(y[1]*y[2] > .(t)), cex = 0.8, adj = 1)
dev.off()

# Equal coordinate marginals with opposite correlations; CopPIT uses an independent forecast reference.
cairo_pdf(file.path(figure_dir, "laajil_marginal_example.pdf"), width = 10.8,
          height = 4, pointsize = 15, family = "sans")
set.seed(20261008)
n <- 10000
sigma_true <- matrix(c(1, .8, .8, 1), 2)
sigma_forecast <- matrix(c(1, -.8, -.8, 1), 2)
y <- mvtnorm::rmvnorm(n, sigma = sigma_true)
forecast_draws <- mvtnorm::rmvnorm(n, sigma = sigma_forecast)
colnames(y) <- colnames(forecast_draws) <- c("Y1", "Y2")
coordinate_pits <- pnorm(y)

set.seed(20261009)
n_reference <- 50000
reference_draws <- mvtnorm::rmvnorm(n_reference, sigma = sigma_forecast)
reference_joint_cdf <- numeric(n_reference)
for (i in seq_len(n_reference)) {
  reference_joint_cdf[i] <- max(0, mvtnorm::pmvnorm(upper = reference_draws[i, ],
    sigma = sigma_forecast, algorithm = mvtnorm::TVPACK()))
}
K_forecast <- ecdf(reference_joint_cdf)
observed_joint_cdf <- numeric(n)
for (i in seq_len(n)) {
  observed_joint_cdf[i] <- max(0, mvtnorm::pmvnorm(upper = y[i, ],
    sigma = sigma_forecast, algorithm = mvtnorm::TVPACK()))
}
coppit <- K_forecast(observed_joint_cdf)

observed_color <- "#167D8D"
forecast_color <- "#C96539"
ink <- "#19324A"
layout(matrix(c(1, 2, 3, 1, 4, 5), 2, byrow = TRUE),
       widths = c(1.12, 1, 1))
par(mar = c(2.5, 2.15, 1.55, .45), mgp = c(1.3, .35, 0),
    tcl = -.2, las = 1, xaxs = "i", yaxs = "i", bty = "l", bg = "white",
    fg = ink, col.axis = ink, col.lab = ink, col.main = ink,
    cex.axis = .78, cex.lab = .83, cex.main = .94)
alpha <- seq(0, 1, length.out = 1001)

par(pty = "s")
plot(NA, xlim = c(-3.8, 3.8), ylim = c(-3.8, 3.8), asp = 1, axes = FALSE,
     xlab = expression(Y[1]), ylab = expression(Y[2]), main = "Samples")
idx <- seq_len(min(1000, n))
points(forecast_draws[idx, ], pch = 16, cex = .42,
       col = adjustcolor(forecast_color, alpha.f = .45))
points(y[idx, ], pch = 16, cex = .42,
       col = adjustcolor(observed_color, alpha.f = .45))
axis(1, at = c(-3, 0, 3))
axis(2, at = c(-3, 0, 3))
box(bty = "l")
legend("topleft", c("Observed", "Forecast"), col = c(observed_color, forecast_color),
       pch = 16, pt.cex = .7, cex = .78, bty = "n", inset = .005,
       y.intersp = .9)
par(pty = "m")

x <- seq(-3.5, 3.5, length.out = 1001)
plot(x, dnorm(x), type = "l", lwd = 2.8, col = observed_color, ylim = c(0, .43),
     xlab = "Coordinate value", ylab = "Density", main = "Marginal Densities",
     axes = FALSE)
lines(x, dnorm(x), col = forecast_color, lwd = 1.8, lty = 2)
axis(1, at = c(-3, 0, 3))
axis(2, at = c(0, .2, .4))
box(bty = "l")

plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = "CDF of marginal PITs",
     xlab = "", ylab = "")
abline(0, 1, col = "#9BA6AF", lty = 3, lwd = 1.6)
axis(1, at = c(0, .5, 1))
axis(2, at = c(0, .5, 1))
box(bty = "l")
lines(alpha, ecdf(coordinate_pits[, 1])(alpha), col = observed_color, lwd = 2)
lines(alpha, ecdf(coordinate_pits[, 2])(alpha), col = ink, lwd = 1.6, lty = 2)
legend("bottomright", expression(Y[1], Y[2]), col = c(observed_color, ink),
       lty = c(1, 2), lwd = c(2, 1.6), bty = "n", cex = .72,
       seg.len = 1.2, y.intersp = .85, inset = .01)

plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
     xlab = expression("Forecast CDF value, " * t), ylab = "Cumulative probability",
     main = "Distribution of forecast CDF values")
lines(alpha, ecdf(observed_joint_cdf)(alpha), col = observed_color, lwd = 2.3)
lines(alpha, K_forecast(alpha), col = forecast_color, lwd = 2.3, lty = 2)
axis(1, at = c(0, .5, 1))
axis(2, at = c(0, .5, 1))
box(bty = "l")
legend("bottomright", c("Observed", "Forecast (Kendall)"),
       col = c(observed_color, forecast_color), lty = c(1, 2), lwd = 2.3,
       bty = "n", cex = .69, seg.len = 1.2, y.intersp = .85,
       inset = .01)

plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = "CDF of CopPIT",
     xlab = "", ylab = "")
abline(0, 1, col = "#9BA6AF", lty = 3, lwd = 1.6)
axis(1, at = c(0, .5, 1))
axis(2, at = c(0, .5, 1))
box(bty = "l")
lines(alpha, ecdf(coppit)(alpha), col = observed_color, lwd = 2.3)
dev.off()
laajil_coppit <- coppit
laajil_y <- y
laajil_reference_joint_cdf <- reference_joint_cdf

# Antisymmetric density on the unit square; swap-invariant scores miss this error.
n <- 12000L
set.seed(20261005)
chunks <- list()
accepted <- 0L
while (accepted < n) {
  batch <- max(1000L, ceiling(2.2 * (n - accepted)))
  proposal <- matrix(runif(2 * batch), ncol = 2)
  keep <- runif(batch) <= (1 + sin(2 * pi * (proposal[, 1] - proposal[, 2]))) / 2
  chunks[[length(chunks) + 1L]] <- proposal[keep, , drop = FALSE]
  accepted <- accepted + sum(keep)
}
y <- do.call(rbind, chunks)[seq_len(n), , drop = FALSE]
difference <- y[, 1] - y[, 2]
average <- rowMeans(y)
z_mean <- ifelse(average <= .5, 2 * average^2, 1 - 2 * (1 - average)^2)
z_spread <- 2 * abs(difference) - difference^2
exceed_count <- rowSums(y > .75)
z_exceed <- pbinom(exceed_count - 1, 2, .25) +
  runif(n) * dbinom(exceed_count, 2, .25)
z_contrast <- ifelse(difference < 0, (difference + 1)^2 / 2,
                     1 - (1 - difference)^2 / 2)

square_y <- y
product <- y[, 1] * y[, 2]
square_coppit <- product * (1 - log(product))
simple_pits <- cbind(Mean = z_mean, Spread = z_spread,
                     Exceedance = z_exceed, Contrast = z_contrast)
cairo_pdf(file.path(figure_dir, "pit_prerank_comparison.pdf"), width = 11,
          height = 2.8, family = "sans", pointsize = 14)
par(mfrow = c(1, 4))
par(mar = c(2.7, 3.2, 1.7, .5), mgp = c(2.0, 0.6, 0), las = 1,
    bg = "white", fg = "#19324A", col.axis = "#19324A",
    col.lab = "#19324A", col.main = "#19324A",
    family = "sans", xaxs = "i", yaxs = "i", bty = "l")
for (j in seq_len(ncol(simple_pits))) {
  hist(simple_pits[, j], breaks = seq(0, 1, length.out = 16),
       freq = FALSE, col = if (j < 4) "#C7E5DF" else "#E7C3B2",
       border = "white", xlim = c(0, 1), ylim = c(0, 2.15),
       main = colnames(simple_pits)[j],
       xlab = expression(Z[g]), ylab = if (j == 1) "Density" else "",
       axes = FALSE)
  abline(h = 1, col = "#19324A", lty = 2, lwd = 2)
  if (j == 4) {
    z <- seq(.001, .999, length.out = 1000)
    contrast_quantile <- ifelse(z <= .5, sqrt(2 * z) - 1, 1 - sqrt(2 * (1 - z)))
    lines(z, 1 + sin(2 * pi * contrast_quantile), col = "#C96539", lwd = 2.5)
  }
  axis(1, at = c(0, .5, 1))
  axis(2, at = 0:2)
  box(bty = "l")
}
dev.off()

# The shared examples compare every diagnostic on the same observations.
set.seed(2042)
z <- matrix(rnorm(4L * n), ncol = 4L)
gaussian_y <- cbind(Y1 = z[, 1], Y2 = .8 * z[, 1] + sqrt(1 - .8^2) * z[, 2])
ua <- pnorm(gaussian_y)
product <- ua[, 1] * ua[, 2]
radial <- pchisq(rowSums(gaussian_y^2), df = 2)
pits_a <- cbind(coordinate_1 = ua[, 1], coordinate_2 = ua[, 2],
  cop_pit = product * (1 - log(product)), mean = pnorm(rowSums(gaussian_y) / sqrt(2)),
  contrast = pnorm((gaussian_y[, 1] - gaussian_y[, 2]) / sqrt(2)),
  ds = radial, hdr = radial, es = radial, latent = radial)
q <- rowSums((square_y - .5)^2)
ds <- pi * q
upper <- q > .25
ds[upper] <- pi * q[upper] - 4 * q[upper] * acos(1 / (2 * sqrt(q[upper]))) +
  2 * sqrt(q[upper] - .25)
set.seed(20261006)
hdr <- runif(n)
pits_b <- cbind(coordinate_1 = square_y[, 1], coordinate_2 = square_y[, 2],
  cop_pit = square_coppit, mean = simple_pits[, "Mean"],
  contrast = simple_pits[, "Contrast"], ds = ds, hdr = hdr,
  latent = pchisq(rowSums(qnorm(square_y)^2), df = 2))
cols <- c(navy = "#17324D", teal = "#238F8B", orange = "#D57732",
          pale_teal = "#C7E5DF", pale_blue = "#DCE5ED", grey = "#626D77")

cairo_pdf(file.path(figure_dir, "comparison_cases.pdf"), width = 8,
          height = 3.2, family = "sans", pointsize = 14)
fill <- colorRampPalette(c("#FFFFFF", "#B3D8D9", cols["teal"]))(81)
par(mfrow = c(1, 2), pty = "s")
par(mar = c(3.05, 2.6, 1.65, .7), mgp = c(1.9, .55, 0), las = 1,
    bg = "white", fg = cols["navy"], col.axis = cols["navy"],
    col.lab = cols["navy"], col.main = cols["navy"],
    cex.axis = .88, cex.lab = .92, cex.main = 1,
    xaxs = "i", yaxs = "i")
grid_a <- seq(-3, 3, length.out = 181)
rho <- .8
xx <- matrix(grid_a, length(grid_a), length(grid_a))
yy <- t(xx)
density_a <- exp(-(xx^2 - 2 * rho * xx * yy + yy^2) / (2 * (1 - rho^2))) /
  (2 * pi * sqrt(1 - rho^2))
image(grid_a, grid_a, density_a, col = fill, useRaster = TRUE,
      xlab = expression(y[1]), ylab = expression(y[2]),
      main = "A: Gaussian dependence error", axes = FALSE, asp = 1)
angle <- seq(0, 2 * pi, length.out = 361)
radius <- sqrt(qchisq(.9, 2))
lines(radius * cos(angle), radius * sin(angle),
      col = cols["orange"], lty = 2, lwd = 2)
axis(1, at = c(-3, 0, 3))
axis(2, at = c(-3, 0, 3))
box(bty = "o")
legend("topleft", "Forecast 90% region", col = cols["orange"],
       lty = 2, lwd = 2, bty = "n", cex = .64, inset = .01)

grid_b <- seq(0, 1, length.out = 181)
density_b <- 1 + sin(2 * pi * outer(grid_b, grid_b, "-"))
image(grid_b, grid_b, density_b, col = fill, zlim = c(0, 2),
      useRaster = TRUE, xlab = expression(y[1]), ylab = expression(y[2]),
      main = "B: Antisymmetric dependence error", axes = FALSE, asp = 1)
axis(1, at = c(0, .5, 1))
axis(2, at = c(0, .5, 1))
box(bty = "o")

dev.off()

for (case in c("gaussian", "square")) {
  pits <- if (case == "gaussian") pits_a else pits_b
  cairo_pdf(file.path(figure_dir, paste0("comparison_", case, "_reliability.pdf")),
            width = 11, height = 3.2, family = "sans", pointsize = 14)
  par(mfrow = c(1, 4))
  par(mar = c(3.05, 3.05, 1.6, .5), mgp = c(1.9, .55, 0), las = 1,
      bg = "white", fg = cols["navy"], col.axis = cols["navy"],
      col.lab = cols["navy"], col.main = cols["navy"],
      cex.axis = .88, cex.lab = .92, cex.main = 1,
      xaxs = "i", yaxs = "i")
  columns <- list(c("coordinate_1", "coordinate_2"), "cop_pit", c("mean", "contrast"),
                  if (case == "gaussian") "ds" else c("ds", "hdr", "latent"))
  labels <- list(c("Coordinate 1", "Coordinate 2"), "CopPIT", c("Mean", "Contrast"),
                 if (case == "gaussian") "HDR / DS / ES / latent" else
                   c("DS", "HDR (jitter)", "Latent radius"))
  titles <- c("Coordinate PITs", "CopPIT", "Simple pre-ranks", "Scores / latent radius")
  grid <- seq(0, 1, length.out = 1001)
  for (panel in seq_along(columns)) {
    values <- pits[, columns[[panel]], drop = FALSE]
    colors <- cols[c("teal", "orange", "navy")][seq_len(ncol(values))]
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         main = titles[panel], xlab = "Predicted proportion",
         ylab = if (panel == 1) "Observed proportion" else "")
    abline(0, 1, col = cols["grey"], lty = 2, lwd = 1.6)
    for (j in seq_len(ncol(values))) {
      lines(grid, ecdf(values[, j])(grid), col = colors[j], lwd = 2,
            lty = if (j == 3L) 3 else 1)
    }
    axis(1, at = c(0, .5, 1))
    axis(2, at = c(0, .5, 1))
    box(bty = "l")
    legend("bottomright", labels[[panel]], col = colors, lwd = 2,
           lty = ifelse(seq_along(labels[[panel]]) == 3L, 3, 1),
           bty = "n", cex = .8, inset = c(.015, .01),
           y.intersp = 1.05, seg.len = 1.4)
  }
  dev.off()
}

# Forecast regions under the opposite-correlation alternative.
alpha <- .80
rho_forecast <- -.8
rho_true <- .8
sigma_forecast <- matrix(c(1, rho_forecast, rho_forecast, 1), 2)
sigma_true <- matrix(c(1, rho_true, rho_true, 1), 2)
set.seed(20261008)
observed <- mvtnorm::rmvnorm(10000, sigma = sigma_true)
forecast_draws <- mvtnorm::rmvnorm(10000, sigma = sigma_forecast)
pc1 <- c(1, -1) / sqrt(2)

# Fixed Gaussian quadrature integrates the forecast-defined boundary.
# Set x = qnorm(cutoff) + s^2 to resolve its infinite left endpoint smoothly.
quadrature <- pracma::gaussLegendre(192, 0, 1)
cutoff_low <- .001
cutoff_high <- .5
for (iteration in seq_len(36L)) {
  cutoff <- (cutoff_low + cutoff_high) / 2
  left <- qnorm(cutoff)
  ss <- sqrt(10 - left) * quadrature$x
  xx <- left + ss^2
  boundary_values <- numeric(length(xx))
  for (j in seq_along(xx)) {
    lo <- left - 1
    hi <- 15
    for (k in seq_len(42L)) {
      mid <- (lo + hi) / 2
      probability <- mvtnorm::pmvnorm(upper = c(xx[j], mid),
        sigma = sigma_forecast, algorithm = mvtnorm::TVPACK())
      if (probability < cutoff) lo <- mid else hi <- mid
    }
    boundary_values[j] <- (lo + hi) / 2
  }
  weights <- quadrature$w * sqrt(10 - left) * 2 * ss * dnorm(xx)
  probability <- cutoff + sum(weights * pnorm(
    (boundary_values - rho_forecast * xx) / sqrt(1 - rho_forecast^2)))
  if (probability < alpha) cutoff_low <- cutoff else cutoff_high <- cutoff
}
cutoffs <- c(qnorm(alpha), sqrt(1 - rho_forecast) * qnorm(alpha), cutoff)
forecast_probability <- c(alpha, alpha, probability)
actual_probability <- c(alpha, pnorm(cutoffs[2] / sqrt(1 - rho_true)),
  cutoff + sum(weights * pnorm((boundary_values - rho_true * xx) /
                                sqrt(1 - rho_true^2))))
joint_scores <- numeric(nrow(observed))
for (i in seq_len(nrow(observed))) {
  joint_scores[i] <- mvtnorm::pmvnorm(upper = observed[i, ],
    sigma = sigma_forecast, algorithm = mvtnorm::TVPACK())
}
scores <- cbind(observed[, 1], drop(observed %*% pc1), joint_scores)
empirical_coverage <- colMeans(sweep(scores, 2, cutoffs, "<="))
results_dir <- file.path(project_dir, "output", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(pre_rank = c("Coordinate", "PCA (PC1)", "Joint CDF / CopPIT"),
  alpha, cutoff = cutoffs, forecast_probability, actual_probability,
  empirical_coverage, n_observed = nrow(observed)),
  file.path(results_dir, "region_miscalibration_coverage.csv"), row.names = FALSE)

ink <- "#19324A"; teal <- "#167D8D"; rust <- "#C96539"
fill <- "#DDEEF0"
limits <- c(-3.6, 3.6)
x <- seq(limits[1], limits[2], length.out = 601)
lo <- limits[1]
hi <- limits[2]
for (iteration in seq_len(42L)) {
  mid <- (lo + hi) / 2
  probability <- mvtnorm::pmvnorm(upper = c(mid, limits[2]),
    sigma = sigma_forecast, algorithm = mvtnorm::TVPACK())
  if (probability < cutoffs[3]) lo <- mid else hi <- mid
}
entry <- (lo + hi) / 2
x <- sort(c(x, entry))
cop_boundary <- rep(Inf, length(x))
for (j in which(pnorm(x) > cutoffs[3])) {
  lo <- qnorm(cutoffs[3]) - 1
  hi <- 15
  for (iteration in seq_len(42L)) {
    mid <- (lo + hi) / 2
    probability <- mvtnorm::pmvnorm(upper = c(x[j], mid),
      sigma = sigma_forecast, algorithm = mvtnorm::TVPACK())
    if (probability < cutoffs[3]) lo <- mid else hi <- mid
  }
  cop_boundary[j] <- (lo + hi) / 2
}
figure_dir <- file.path(project_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(figure_dir, "region_miscalibration.pdf"), width = 11, height = 3.3,
          pointsize = 16, family = "sans", bg = "white")
par(mfrow = c(1, 3), mar = c(2.3, 2.5, 1.4, 2.5), oma = c(0, 0, 0, 0),
    mgp = c(1.45, .35, 0), tcl = -.17, las = 1, pty = "s",
    xaxs = "i", yaxs = "i", fg = ink, col.axis = ink,
    col.lab = ink, col.main = ink, cex = 1, cex.axis = .86, cex.lab = .95)
titles <- c("Marginal", "PCA (PC1)", "Joint CDF (CopPIT)")
for (i in 1:3) {
  plot(NA, xlim = limits, ylim = limits, asp = 1, axes = FALSE,
       xlab = expression(y[1]), ylab = expression(y[2]))
  if (i == 1) rect(limits[1], limits[1], cutoffs[1], limits[2],
                   col = fill, border = NA)
  if (i == 2) {
    lower <- pmax(limits[1], pmin(limits[2], x - sqrt(2) * cutoffs[2]))
    polygon(c(x, rev(x)), c(lower, rep(limits[2], length(x))),
            col = fill, border = NA)
  }
  if (i == 3) {
    upper <- pmax(limits[1], pmin(limits[2], cop_boundary))
    polygon(c(x, rev(x)), c(upper, rep(limits[1], length(x))),
            col = fill, border = NA)
  }
  abline(h = 0, v = 0, col = "#DDE1E4", lwd = .7)
  points(forecast_draws[1:500, ], pch = 16, cex = .43,
         col = adjustcolor(rust, alpha.f = .45))
  points(observed[1:500, ], pch = 16, cex = .43,
         col = adjustcolor(teal, alpha.f = .45))
  if (i == 1) abline(v = cutoffs[1], col = teal, lwd = 1.5)
  if (i == 2) abline(a = -sqrt(2) * cutoffs[2], b = 1, col = teal, lwd = 1.5)
  if (i == 3) {
    visible <- x > entry
    lines(c(entry, x[visible]), c(limits[2], cop_boundary[visible]),
          col = teal, lwd = 1.5)
  }
  axis(1, at = c(-3, 0, 3)); axis(2, at = c(-3, 0, 3))
  box(col = "#7B8A96", lwd = .6)
  mtext(titles[i], side = 3, line = .35, font = 2, cex = .92)
}
dev.off()

# 3. Forecast samples and pre-rank regions ---------------------------------

message("Pre-rank regions")
grid_axis <- seq(-4.5, 4.5, length.out = 301L)
grid <- as.matrix(expand.grid(y1 = grid_axis, y2 = grid_axis))
u_axis <- seq(0, 1, length.out = 301L)
u_grid <- as.matrix(expand.grid(u1 = u_axis, u2 = u_axis))
z_axis <- seq(-3.5, 3.5, length.out = 301L)
z_grid <- as.matrix(expand.grid(z1 = z_axis, z2 = z_axis))
quadrature <- pracma::gaussLegendre(384L, a = -1, b = 1)
forecasts <- vector("list", 3L)

for (row in 1:3) {
  set.seed(42)
  sigma <- .65
  covariance <- switch(row, diag(c(4 + sigma^2, sigma^2)),
    matrix(c(1, -.8, -.8, 1), 2), diag(c(1, 1.44)))
  draws <- vector("list", 2L)
  for (sample_index in 1:2) {
    n <- c(1000L, 5000L)[sample_index]
    if (row == 1L) {
      draws[[sample_index]] <- cbind(
        sample(c(-2, 2), n, replace = TRUE) + rnorm(n, sd = sigma),
        rnorm(n, sd = sigma))
    } else {
      z <- matrix(rnorm(2L * n), ncol = 2L)
      draws[[sample_index]] <- if (row == 2L) z %*% chol(covariance) else
        cbind(z[, 1], .8 * (z[, 1]^2 - 1) + .4 * z[, 2])
    }
  }
  reference <- draws[[1L]]
  base_draws <- draws[[2L]]
  y <- rbind(base_draws, grid)
  eig <- eigen(covariance, symmetric = TRUE)
  direction <- eig$vectors[, 1L]
  pivot <- which(abs(direction) >= max(abs(direction)) - 1e-12)[1L]
  if (direction[pivot] < 0) direction <- -direction

  if (row == 1L) {
    marginal_pits <- cbind((pnorm(y[, 1], -2, sigma) +
      pnorm(y[, 1], 2, sigma)) / 2, pnorm(y[, 2], 0, sigma))
    joint_cdf <- marginal_pits[, 1] * marginal_pits[, 2]
    density <- (dnorm(y[, 1], -2, sigma) + dnorm(y[, 1], 2, sigma)) / 2 *
      dnorm(y[, 2], 0, sigma)
    latent <- qnorm(pmin(1 - 1e-15, pmax(1e-15, marginal_pits)))
    dim(latent) <- dim(marginal_pits)
    pit_joint_cdf <- u_grid[, 1] * u_grid[, 2]
  } else if (row == 2L) {
    marginal_pits <- pnorm(y)
    density <- mvtnorm::dmvnorm(y, sigma = covariance)
    whitening <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
    latent <- y %*% whitening
    interior <- pmin(1 - 1e-12, pmax(1e-12, u_axis))
    pit_coordinates <- as.matrix(expand.grid(qnorm(interior), qnorm(interior)))
    cdf_points <- rbind(y, pit_coordinates)
    cdf_values <- numeric(nrow(cdf_points))
    for (i in seq_len(nrow(cdf_points))) {
      cdf_values[i] <- max(0, mvtnorm::pmvnorm(upper = cdf_points[i, ],
        sigma = covariance, algorithm = mvtnorm::TVPACK()))
    }
    joint_cdf <- cdf_values[seq_len(nrow(y))]
    pit_joint_cdf <- cdf_values[nrow(y) + seq_len(nrow(u_grid))]
  } else {
    density <- dnorm(y[, 1]) * dnorm(y[, 2], .8 * (y[, 1]^2 - 1), .4)
    latent <- cbind(y[, 1], (y[, 2] - .8 * (y[, 1]^2 - 1)) / .4)

    # Invert the banana's second marginal CDF at the PIT-grid probabilities.
    probabilities <- pmin(1 - 1e-12, pmax(1e-12, u_axis))
    lower <- rep(-16, length(probabilities))
    upper <- rep(64, length(probabilities))
    nodes <- 9 * quadrature$x
    mass <- 9 * quadrature$w * dnorm(nodes)
    for (iteration in 1:45) {
      middle <- (lower + upper) / 2
      cdf <- drop(pnorm(outer(middle, .8 * (nodes^2 - 1), "-") / .4) %*% mass)
      below <- cdf < probabilities
      lower[below] <- middle[below]
      upper[!below] <- middle[!below]
    }
    pit_coordinates <- as.matrix(expand.grid(qnorm(probabilities), (lower + upper) / 2))

    # Integrate phi(z) Phi((y2 - .8*(z^2-1))/.4) over z <= y1.
    # The [-9,9] truncation omits less than 2.3e-19 Gaussian probability.
    marginal_y2 <- unique(y[, 2])
    cdf_points <- rbind(y, pit_coordinates, cbind(Inf, marginal_y2))
    cdf_values <- numeric(nrow(cdf_points))
    for (start in seq.int(1L, nrow(cdf_points), by = 2048L)) {
      index <- start:min(start + 2047L, nrow(cdf_points))
      half_width <- (pmin(cdf_points[index, 1], 9) + 9) / 2
      nodes <- outer(half_width, quadrature$x + 1) - 9
      integrand <- dnorm(nodes) *
        pnorm((cdf_points[index, 2] - .8 * (nodes^2 - 1)) / .4)
      cdf_values[index] <- pmin(pnorm(cdf_points[index, 1]),
        pmax(0, half_width * drop(integrand %*% quadrature$w)))
    }
    joint_cdf <- cdf_values[seq_len(nrow(y))]
    pit_joint_cdf <- cdf_values[nrow(y) + seq_len(nrow(u_grid))]
    marginal_cdf <- tail(cdf_values, length(marginal_y2))
    marginal_pits <- cbind(pnorm(y[, 1]), marginal_cdf[match(y[, 2], marginal_y2)])
  }
  pit_joint_cdf[u_grid[, 1] == 0 | u_grid[, 2] == 0] <- 0
  pit_joint_cdf[u_grid[, 1] == 1] <- u_grid[u_grid[, 1] == 1, 2]
  pit_joint_cdf[u_grid[, 2] == 1] <- u_grid[u_grid[, 2] == 1, 1]

  # The frozen reference defines the energy/variogram scores; separate draws set cutoffs.
  energy <- numeric(nrow(y))
  for (start in seq.int(1L, nrow(y), by = 512L)) {
    index <- start:min(start + 511L, nrow(y))
    energy[index] <- rowMeans(sqrt(outer(y[index, 1], reference[, 1], "-")^2 +
      outer(y[index, 2], reference[, 2], "-")^2))
  }
  standardized <- sweep(y, 2L, sqrt(diag(covariance)), "/")
  scores <- cbind(
    coordinate = y[, 1], location = rowMeans(y), pca_projection = drop(y %*% direction),
    radius = sqrt(rowSums(y^2)),
    mahalanobis = sqrt(pmax(0, mahalanobis(y, c(0, 0), covariance))),
    box = pmax(abs(standardized[, 1]), abs(standardized[, 2])),
    orthant = pmax(standardized[, 1], standardized[, 2]),
    joint_cdf = joint_cdf, average_rank = rowMeans(marginal_pits),
    band_depth = -rowMeans(marginal_pits * (1 - marginal_pits)),
    hdr = density, energy = energy, scale = rowMeans((y - rowMeans(y))^2),
    variogram = (abs(y[, 1] - y[, 2]) - mean(abs(reference[, 1] - reference[, 2])))^2,
    latent = sqrt(rowSums(latent^2)))
  forecasts[[row]] <- list(base_draws = base_draws, direction = direction,
    sample_scores = scores[1:5000, ], grid_scores = scores[-(1:5000), ],
    marginal_pits_sample = marginal_pits[1:5000, ], latent_sample = latent[1:5000, ],
    pit_scores = cbind(joint_cdf = pit_joint_cdf,
      average_rank = rowMeans(u_grid), band_depth = -rowMeans(u_grid * (1 - u_grid))))
  stopifnot(all(is.finite(scores)), all(joint_cdf >= 0 & joint_cdf <= 1))
}

# Render only the region views used by the main deck and its supplements.
region_methods <- c("coordinate", "location", "pca_projection", "radius", "mahalanobis",
  "box", "orthant", "joint_cdf", "average_rank", "band_depth", "hdr", "energy",
  "scale", "variogram")
plots <- data.frame(
  name = c(paste0(region_methods, "_alpha80"), "hdr_alpha50", "hdr_alpha95",
    "pit_joint_cdf_alpha80", "pit_average_rank_alpha80", "pit_band_depth_alpha80",
    "latent_space_alpha80", "pca_gaussian_levels"),
  method = c(region_methods, "hdr", "hdr", "joint_cdf", "average_rank", "band_depth",
    "latent", "pca_projection"),
  coordinates = c(rep("original", 16), rep("pit", 3), "latent", "original"),
  alpha = c(rep(.8, 14), .5, .95, rep(.8, 5)))
ink <- "#19324A"
teal <- "#167D8D"
rust <- "#C96539"
fill <- "#DDEEF0"
dots <- adjustcolor(ink, alpha.f = .37)

for (plot_index in seq_len(nrow(plots))) {
  method <- plots$method[plot_index]
  coordinates <- plots$coordinates[plot_index]
  pca_levels <- plots$name[plot_index] == "pca_gaussian_levels"
  cairo_pdf(file.path(figure_dir, paste0("prerank_", plots$name[plot_index], ".pdf")),
    width = 11, height = if (pca_levels) 6.1 else 3.3,
    pointsize = if (pca_levels) 14 else 16, family = "sans", bg = "white")
  if (pca_levels) {
    par(mfrow = c(3, 3), mar = c(1.5, 1.45, 1.35, .2), oma = c(0, 1.5, 0, 0),
      mgp = c(.90, .25, 0), tcl = -.13, cex.axis = .72, cex.lab = .8)
  } else {
    par(mfrow = c(1, 3), mar = c(2.3, 2.5, 1.4, 2.5), oma = rep(0, 4),
      mgp = c(1.45, .35, 0), tcl = -.17, cex.axis = .86, cex.lab = .95)
  }
  par(las = 1, pty = "s", xaxs = "i", yaxs = "i", cex = 1,
    bg = "white", fg = ink, col.axis = ink, col.lab = ink, col.main = ink)
  for (row in 1:3) {
    f <- forecasts[[row]]
    levels <- if (pca_levels) c(.5, .8, .95) else plots$alpha[plot_index]
    for (column in seq_along(levels)) {
      alpha <- levels[column]
      cutoff <- unname(quantile(f$sample_scores[, method],
        if (method == "hdr") 1 - alpha else alpha))
      if (coordinates == "original") {
        axis_values <- grid_axis
        values <- matrix(f$grid_scores[, method], nrow = length(axis_values))
        point_data <- f$base_draws[1:300, ]
        ticks <- c(-4, 0, 4)
        xlab <- expression(y[1]); ylab <- expression(y[2])
      } else if (coordinates == "pit") {
        axis_values <- u_axis
        values <- matrix(f$pit_scores[, method], nrow = length(axis_values))
        point_data <- f$marginal_pits_sample[1:300, ]
        ticks <- c(0, .5, 1)
        xlab <- expression(u[1]); ylab <- expression(u[2])
      } else {
        axis_values <- z_axis
        values <- matrix(sqrt(rowSums(z_grid^2)), nrow = length(axis_values))
        point_data <- f$latent_sample[1:300, ]
        ticks <- c(-3, 0, 3)
        xlab <- expression(z[1]); ylab <- expression(z[2])
      }
      if (pca_levels && row != 3L) xlab <- ""
      if (pca_levels && column != 1L) ylab <- ""
      inside <- if (method == "hdr") values >= cutoff else values <= cutoff
      plot(NA, xlim = range(axis_values), ylim = range(axis_values), asp = 1,
        axes = FALSE, xlab = xlab, ylab = ylab)
      # Consecutive horizontal runs keep the shaded PDFs fully vector based.
      half_cell <- diff(axis_values)[1L] / 2
      for (j in seq_along(axis_values)) {
        runs <- rle(inside[, j])
        last <- cumsum(runs$lengths)
        first <- last - runs$lengths + 1L
        selected <- which(runs$values)
        if (length(selected)) rect(axis_values[first[selected]] - half_cell, axis_values[j] - half_cell,
          axis_values[last[selected]] + half_cell, axis_values[j] + half_cell,
          col = fill, border = NA)
      }
      abline(h = 0, v = 0, col = "#DAE0E5", lwd = .55)
      points(point_data, pch = 16, cex = .42, col = dots)
      contour(axis_values, axis_values, values, levels = cutoff, add = TRUE,
        drawlabels = FALSE, col = teal, lwd = 1.5)
      axis(1, at = ticks); axis(2, at = ticks)
      box(col = "#7B8A96", lwd = .6)
      title <- if (pca_levels) sprintf("PC1: %.0f%% region", 100 * alpha) else LETTERS[row]
      mtext(title, side = 3, line = .35, font = 2, cex = .92)
      if (method == "pca_projection") arrows(0, 0, 2.7 * f$direction[1],
        2.7 * f$direction[2], length = .075, lwd = 1.4, col = rust)
    }
    if (pca_levels) mtext(LETTERS[row], side = 2, outer = TRUE,
      at = 1 - (row - .5) / 3, line = .1, las = 1, font = 2, cex = 1)
  }
  dev.off()
}

# Normalized dependence needs d=3: in d=2 the defined score is constant.
slice_axis <- seq(-3, 3, length.out = 301L)
xy <- as.matrix(expand.grid(y1 = slice_axis, y2 = slice_axis))
dependence_cutoffs <- numeric(3)
for (row in 1:3) {
  set.seed(42)
  z <- matrix(rnorm(30000L), ncol = 3L)
  covariance <- switch(row, diag(3), .8^abs(outer(1:3, 1:3, "-")), diag(c(1, 1.44, 1)))
  y <- if (row == 3L) cbind(z[, 1], .8 * (z[, 1]^2 - 1) + .4 * z[, 2], z[, 3]) else
    z %*% chol(covariance)
  variance <- rowMeans((y - rowMeans(y))^2)
  score <- -rowSums((y[, 1:2] - y[, 2:3])^2) / (4 * variance)
  dependence_cutoffs[row] <- unname(quantile(score, .8))
}
for (third in c(0, 1)) {
  y <- cbind(xy, third)
  variance <- rowMeans((y - rowMeans(y))^2)
  score <- -rowSums((y[, 1:2] - y[, 2:3])^2) / (4 * variance)
  score[variance == 0] <- NA_real_
  values <- matrix(score, nrow = length(slice_axis))
  cairo_pdf(file.path(figure_dir, paste0("prerank_dependence_slice", third, "_alpha80.pdf")),
    width = 11, height = 3.3, pointsize = 16, family = "sans", bg = "white")
  par(mfrow = c(1, 3), mar = c(2.3, 2.5, 1.4, 2.5), oma = rep(0, 4),
    mgp = c(1.45, .35, 0), tcl = -.17, las = 1, pty = "s",
    xaxs = "i", yaxs = "i", cex = 1, cex.axis = .86, cex.lab = .95,
    bg = "white", fg = ink, col.axis = ink, col.lab = ink, col.main = ink)
  for (row in 1:3) {
    cutoff <- dependence_cutoffs[row]
    inside <- !is.na(values) & values <= cutoff
    plot(NA, xlim = range(slice_axis), ylim = range(slice_axis), asp = 1,
      axes = FALSE, xlab = expression(y[1]), ylab = expression(y[2]))
    half_cell <- diff(slice_axis)[1L] / 2
    for (j in seq_along(slice_axis)) {
      runs <- rle(inside[, j])
      last <- cumsum(runs$lengths)
      first <- last - runs$lengths + 1L
      selected <- which(runs$values)
      if (length(selected)) rect(slice_axis[first[selected]] - half_cell, slice_axis[j] - half_cell,
        slice_axis[last[selected]] + half_cell, slice_axis[j] + half_cell,
        col = fill, border = NA)
    }
    abline(h = 0, v = 0, col = "#DAE0E5", lwd = .55)
    contour(slice_axis, slice_axis, values, levels = cutoff, add = TRUE,
      drawlabels = FALSE, col = teal, lwd = 1.5)
    axis(1, at = c(-3, 0, 3)); axis(2, at = c(-3, 0, 3))
    box(col = "#7B8A96", lwd = .6)
    mtext(LETTERS[row], side = 3, line = .35, font = 2, cex = .92)
    points(third, third, pch = 4, lwd = 1.5, cex = .8, col = rust)
  }
  dev.off()
}

# 4. Band-depth regions at alpha = 0.2 ------------------------------------

message("Band-depth regions at alpha = 0.2")

alpha <- .2
ink <- "#19324A"
teal <- "#167D8D"
fill <- "#DDEEF0"
dots <- adjustcolor(ink, alpha.f = .37)
labels <- c("A   Bimodal", "B   Correlated Gaussian", "C   Banana")
quadrature <- pracma::gaussLegendre(384L, a = -1, b = 1)
banana_nodes <- 9 * quadrature$x
banana_weights <- 9 * quadrature$w * dnorm(banana_nodes)
theta <- seq(0, 2 * pi, length.out = 1601L)
band_boundaries_y <- band_boundaries_u <- band_shown_y <- band_shown_u <- vector("list", 3L)
band_checks <- vector("list", 3L)

for (row in 1:3) {
  # Consume the original 1,000 reference draws before the 5,000 cutoff draws.
  set.seed(42L)
  for (sample_n in c(1000L, 5000L)) {
    if (row == 1L) {
      y <- cbind(sample(c(-2, 2), sample_n, replace = TRUE) + rnorm(sample_n, sd = .65),
                 rnorm(sample_n, sd = .65))
    } else {
      z <- matrix(rnorm(2L * sample_n), ncol = 2L)
      if (row == 2L) y <- z %*% chol(matrix(c(1, -.8, -.8, 1), 2L))
      if (row == 3L) y <- cbind(z[, 1], .8 * (z[, 1]^2 - 1) + .4 * z[, 2])
    }
  }
  cutoff_draws <- y
  set.seed(7020L + row)
  if (row == 1L) {
    fresh <- cbind(sample(c(-2, 2), 20000L, replace = TRUE) + rnorm(20000L, sd = .65),
                   rnorm(20000L, sd = .65))
  } else {
    z <- matrix(rnorm(40000L), ncol = 2L)
    if (row == 2L) fresh <- z %*% chol(matrix(c(1, -.8, -.8, 1), 2L))
    if (row == 3L) fresh <- cbind(z[, 1], .8 * (z[, 1]^2 - 1) + .4 * z[, 2])
  }
  y <- rbind(cutoff_draws, fresh)
  if (row == 1L) u <- cbind((pnorm(y[, 1], -2, .65) + pnorm(y[, 1], 2, .65))/2,
                            pnorm(y[, 2], 0, .65))
  if (row == 2L) u <- pnorm(y)
  if (row == 3L) {
    u <- cbind(pnorm(y[, 1]), numeric(nrow(y)))
    for (start in seq.int(1L, nrow(y), by = 2048L)) {
      index <- start:min(start + 2047L, nrow(y))
      integrand <- pnorm(outer(y[index, 2], .8 * (banana_nodes^2 - 1), "-") / .4)
      u[index, 2] <- pmin(1, pmax(0, drop(integrand %*% banana_weights)))
    }
  }
  scores <- -rowMeans(u * (1 - u))
  cutoff <- unname(quantile(scores[1:5000], alpha))
  radius <- sqrt(2 * (cutoff + .25))
  pit_scores <- rowSums((u[1:5000, ] - .5)^2) / 2 - .25
  band_boundaries_u[[row]] <- cbind(.5 + radius * cos(theta), .5 + radius * sin(theta))
  boundary <- band_boundaries_u[[row]]
  original <- qnorm(boundary)
  if (row == 1L) original[, 2] <- qnorm(boundary[, 2], 0, .65)
  if (row != 2L) {
    coordinate <- if (row == 1L) 1L else 2L
    probability <- boundary[, coordinate]
    lower <- rep(-10, length(probability))
    upper <- rep(10, length(probability))
    # Vectorized bisection inverts the mixture and banana marginal CDFs.
    for (iteration in 1:40) {
      midpoint <- (lower + upper) / 2
      if (row == 1L) values <- (pnorm(midpoint, -2, .65) + pnorm(midpoint, 2, .65))/2
      if (row == 3L) {
        integrand <- pnorm(outer(midpoint, .8 * (banana_nodes^2 - 1), "-") / .4)
        values <- pmin(1, pmax(0, drop(integrand %*% banana_weights)))
      }
      below <- values < probability
      lower[below] <- midpoint[below]
      upper[!below] <- midpoint[!below]
    }
    original[, coordinate] <- (lower + upper) / 2
  }
  band_boundaries_y[[row]] <- original
  band_shown_y[[row]] <- cutoff_draws[1:300, ]
  band_shown_u[[row]] <- u[1:300, ]
  band_checks[[row]] <- data.frame(
    forecast = LETTERS[row], alpha = alpha, cutoff = cutoff, pit_radius = radius,
    cutoff_sample_coverage = mean(scores[1:5000] <= cutoff),
    independent_coverage = mean(scores[5001:25000] <= cutoff),
    max_score_difference = max(abs(scores[1:5000] - pit_scores)),
    membership_mismatches = sum((scores[1:5000] <= cutoff) !=
                                (rowSums((u[1:5000, ] - .5)^2) <= radius^2)),
    all_original_coordinate_bounds_finite = all(is.finite(original)))
}
band_checks <- do.call(rbind, band_checks)
stopifnot(all(band_checks$membership_mismatches == 0L),
          all(band_checks$pit_radius < .5),
          all(band_checks$all_original_coordinate_bounds_finite),
          all(abs(band_checks$independent_coverage - alpha) < .015))

cairo_pdf(file.path(figure_dir, "prerank_band_depth_alpha20_six.pdf"),
          width = 11, height = 5, pointsize = 16, family = "sans", bg = "white")
par(oma = rep(0, 4), bg = "white", fg = ink, col.axis = ink, col.lab = ink, col.main = ink)
layout(rbind(rep(1L, 3L), 2:4, rep(5L, 3L), 6:8), heights = c(.13, 1, .13, 1))
for (pit in c(FALSE, TRUE)) {
  par(mar = rep(0, 4), pty = "m", cex = 1, xaxs = "i", yaxs = "i")
  plot.new()
  text(.028, .55, if (pit) "Marginal PIT coordinates" else "Original coordinates",
       adj = c(0, .5), font = 2, cex = .9, col = ink)
  for (row in 1:3) {
    boundary <- if (pit) band_boundaries_u[[row]] else band_boundaries_y[[row]]
    point_data <- if (pit) band_shown_u[[row]] else band_shown_y[[row]]
    par(mar = c(1.65, 2.25, .8, 2.25), mgp = c(.95, .22, 0), tcl = -.17,
        las = 1, pty = "s", xaxs = "i", yaxs = "i", cex = 1, cex.axis = .8, cex.lab = .86)
    limits <- if (pit) c(0, 1) else c(-4.5, 4.5)
    plot(NA, xlim = limits, ylim = limits, asp = 1, axes = FALSE,
         xlab = if (pit) expression(u[1]) else expression(y[1]),
         ylab = if (pit) expression(u[2]) else expression(y[2]))
    polygon(boundary, col = fill, border = NA)
    abline(h = if (pit) .5 else 0, v = if (pit) .5 else 0, col = "#DAE0E5", lwd = .55)
    points(point_data, pch = 16, cex = .36, col = dots)
    lines(boundary, col = teal, lwd = 1.5)
    ticks <- if (pit) c(0, .5, 1) else c(-4, 0, 4)
    axis(1, at = ticks)
    axis(2, at = ticks)
    box(col = "#7B8A96", lwd = .6)
    if (!pit) mtext(labels[row], side = 3, line = .25, font = 2, cex = .88)
  }
}
dev.off()
write.csv(band_checks, file.path(result_dir, "band_depth_alpha20_checks.csv"), row.names = FALSE)

# 5. Population PIT comparison -------------------------------------------

set.seed(20261015L)
n <- 10000L
n_reference <- 50000L
n_null <- 10000L
d <- 10L
distance <- abs(outer(seq_len(d), seq_len(d), "-"))
observations <- matrix(rnorm(n * d), n, d) %*% chol(exp(-distance))
reference_normals <- matrix(rnorm(n_reference * d), n_reference, d)
null_normals <- matrix(rnorm(n_null * d), n_null, d)
scenarios <- data.frame(
  scenario = c("mean_low", "mean_high", "variance_low", "variance_high",
               "correlation_low", "correlation_high"),
  mu = c(-.5, .5, 0, 0, 0, 0), variance = c(1, 1, .85, 1.25, 1, 1),
  tau = c(1, 1, 1, 1, .5, 2))
methods <- c("coppit", "average", "negative_band_depth", "energy",
             "location", "scale", "dependence", "pca", "hdr")
pit_values <- matrix(NA_real_, n * nrow(scenarios), length(methods),
                     dimnames = list(NULL, methods))
reference_tolerance <- sqrt(log(108 / .001) / (2 * n_reference))
null_tolerance <- reference_tolerance + sqrt(log(108 / .001) / (2 * n_null))
checks <- list()
null_tables <- list()

# Standardization preserves score orderings, so reference draws are shared
# by forecasts with the same correlation. Reference and evaluation are independent.
for (tau in c(.5, 1, 2)) {
  message("Population PITs: tau = ", tau)
  sigma <- exp(-distance / tau)
  decomposition <- eigen(sigma, symmetric = TRUE)
  direction <- decomposition$vectors[, 1L]
  pivot <- which(abs(direction) >= max(abs(direction)) - 1e-12)[1L]
  if (direction[pivot] < 0) direction <- -direction
  matching <- which(scenarios$tau == tau)
  z <- rbind(reference_normals %*% chol(sigma), null_normals %*% chol(sigma))
  probe_indices <- seq_len(128L)
  for (i in matching) {
    probe_indices <- c(probe_indices, nrow(z) + seq_len(128L))
    z <- rbind(z, (observations - scenarios$mu[i]) / sqrt(scenarios$variance[i]))
  }
  u <- pnorm(z)
  location <- rowMeans(z)
  scale <- rowMeans((z - location)^2)
  scores <- cbind(
    coppit = NA_real_, average = rowMeans(u),
    negative_band_depth = -rowMeans(u * (1 - u)), energy = NA_real_,
    location = location, scale = scale,
    dependence = -(rowSums((z[, -d] - z[, -1L])^2) / (2 * (d - 1))) / scale,
    pca = drop(z %*% direction), hdr = .5 * mahalanobis(z, rep(0, d), sigma))

  # E||X-z|| from the Gaussian Laplace transform, using 64 quadrature nodes.
  # With t=(u/(a*(1-u)))^2 and a=sqrt(tr(Sigma)), integrate
  # a/sqrt(pi) * (1-L(t))/u^2 on [0,1]. Check 128 nodes on representative points.
  for (energy_nodes in c(64L, 128L)) {
    k <- seq_len(energy_nodes - 1L)
    jacobi <- matrix(0, energy_nodes, energy_nodes)
    jacobi[cbind(k, k + 1L)] <- k / sqrt(4 * k^2 - 1)
    jacobi[cbind(k + 1L, k)] <- k / sqrt(4 * k^2 - 1)
    quadrature <- eigen(jacobi, symmetric = TRUE)
    nodes <- (quadrature$values + 1) / 2
    energy_scale <- sqrt(sum(decomposition$values))
    t_nodes <- (nodes / (energy_scale * (1 - nodes)))^2
    denominators <- 1 + outer(2 * decomposition$values, t_nodes)
    coefficients <- matrix(rep(t_nodes, each = d), d) / denominators
    log_determinant <- -.5 * colSums(log1p(outer(2 * decomposition$values, t_nodes)))
    weights <- energy_scale / sqrt(pi) * quadrature$vectors[1L, ]^2 / nodes^2
    indices <- if (energy_nodes == 64L) seq_len(nrow(z)) else probe_indices
    energy <- numeric(length(indices))
    for (start in seq.int(1L, length(indices), by = 10000L)) {
      rows <- start:min(length(indices), start + 9999L)
      coordinates <- z[indices[rows], , drop = FALSE] %*% decomposition$vectors
      log_laplace <- sweep(-(coordinates^2 %*% coefficients), 2L, log_determinant, "+")
      energy[rows] <- drop((-expm1(log_laplace)) %*% weights)
    }
    if (energy_nodes == 64L) scores[, "energy"] <- energy
    else checks[[length(checks) + 1L]] <- data.frame(
      check = paste0("energy_convergence_tau", tau),
      error = max(abs(energy - scores[probe_indices, "energy"])), tolerance = 1e-9)
  }

  # Deterministic AR(1) joint CDF: integrate the first d-1 coordinates and
  # use the analytic conditional CDF for the last. Log probabilities preserve
  # tail orderings. Batches bound memory; 64 nodes check the 48-node calculation.
  rho <- exp(-1 / tau)
  conditional_sd <- sqrt(1 - rho^2)
  for (cdf_nodes in c(48L, 64L)) {
    quadrature <- pracma::gaussLegendre(cdf_nodes, -1, 1)
    source_index <- rep(seq_len(cdf_nodes), each = cdf_nodes)
    target_index <- rep(seq_len(cdf_nodes), times = cdf_nodes)
    indices <- if (cdf_nodes == 48L) seq_len(nrow(z)) else probe_indices
    log_cdf <- numeric(length(indices))
    for (start in seq.int(1L, length(indices), by = 256L)) {
      rows <- start:min(length(indices), start + 255L)
      upper <- z[indices[rows], , drop = FALSE]
      nb <- nrow(upper)
      lower <- rep(-12, nb)
      for (j in seq_len(d)) lower <- pmin(lower, upper[, j] - 6)
      half <- (pmin(upper[, 1L], 12) - lower) / 2
      x <- outer(half, quadrature$x) + lower + half
      mass <- sweep(dnorm(x) * half, 2L, quadrature$w, "*")
      total <- rowSums(mass)
      log_total <- log(total)
      mass <- mass / total
      for (j in 2:(d - 1L)) {
        half <- (pmin(upper[, j], 12) - lower) / 2
        next_x <- outer(half, quadrature$x) + lower + half
        density <- dnorm((next_x[, target_index] - rho * x[, source_index]) /
                           conditional_sd) / conditional_sd
        density <- rowSums(array(density * mass[, source_index],
                                 c(nb, cdf_nodes, cdf_nodes)), dims = 2L)
        mass <- sweep(density * half, 2L, quadrature$w, "*")
        total <- rowSums(mass)
        log_total <- log_total + log(total)
        mass <- mass / total
        x <- next_x
      }
      log_cdf[rows] <- log_total + log(rowSums(mass *
        pnorm((upper[, d] - rho * x) / conditional_sd)))
    }
    stopifnot(all(is.finite(log_cdf)), max(log_cdf) < 1e-8)
    log_cdf <- pmin(0, log_cdf)
    if (cdf_nodes == 48L) scores[, "coppit"] <- log_cdf
    else checks[[length(checks) + 1L]] <- data.frame(
      check = paste0("joint_cdf_log_convergence_tau", tau),
      error = max(abs(log_cdf - scores[probe_indices, "coppit"])), tolerance = 1e-5)
  }
  stopifnot(all(is.finite(scores)))

  # The empirical reference CDF is the fraction of reference scores <= a score.
  all_pits <- matrix(NA_real_, nrow(z), length(methods), dimnames = list(NULL, methods))
  for (method in methods) {
    reference <- sort(scores[seq_len(n_reference), method])
    all_pits[, method] <- findInterval(scores[, method], reference) / n_reference
    null <- sort(all_pits[n_reference + seq_len(n_null), method])
    checks[[length(checks) + 1L]] <- data.frame(
      check = paste0("null_tau", tau, "_", method),
      error = max(seq_len(n_null) / n_null - null,
                  null - (seq_len(n_null) - 1L) / n_null), tolerance = null_tolerance)
  }
  reference <- sort(scores[seq_len(n_reference), "coppit"])
  checks[[length(checks) + 1L]] <- data.frame(
    check = paste0("coppit_convergence_tau", tau),
    error = max(abs(findInterval(log_cdf, reference) / n_reference -
                    all_pits[probe_indices, "coppit"])), tolerance = 2 / n_reference)
  depth_rows <- n_reference + seq_len(128L)
  checks[[length(checks) + 1L]] <- data.frame(
    check = paste0("negative_depth_orientation_tau", tau),
    error = max(abs(all_pits[depth_rows, "negative_band_depth"] +
      findInterval(-scores[depth_rows, "negative_band_depth"],
                   sort(-scores[seq_len(n_reference), "negative_band_depth"])) /
      n_reference - 1)), tolerance = 1e-12)
  shifted <- 2 * z[probe_indices, ] + 3
  shifted_score <- -(rowSums((shifted[, -d] - shifted[, -1L])^2) /
                     (2 * (d - 1))) / rowMeans((shifted - rowMeans(shifted))^2)
  checks[[length(checks) + 1L]] <- data.frame(
    check = paste0("dependence_shift_scale_invariance_tau", tau),
    error = max(abs(scores[probe_indices, "dependence"] - shifted_score)), tolerance = 1e-12)
  null_tables[[as.character(tau)]] <- data.frame(
    tau = tau, case = seq_len(n_null), all_pits[n_reference + seq_len(n_null), ])
  evaluation_rows <- (n_reference + 1L):nrow(z)
  exact <- cbind(
    location = pnorm(scores[evaluation_rows, "location"], sd = sqrt(mean(sigma))),
    pca = pnorm(scores[evaluation_rows, "pca"],
                sd = sqrt(drop(crossprod(direction, sigma %*% direction)))),
    hdr = pchisq(2 * scores[evaluation_rows, "hdr"], df = d))
  for (method in colnames(exact)) checks[[length(checks) + 1L]] <- data.frame(
    check = paste0("reference_tau", tau, "_", method),
    error = max(abs(all_pits[n_reference + seq_len(n_null), method] -
                    exact[seq_len(n_null), method])), tolerance = reference_tolerance)
  for (j in seq_along(matching)) {
    from <- n_reference + n_null + (j - 1L) * n + seq_len(n)
    to <- (matching[j] - 1L) * n + seq_len(n)
    pit_values[to, ] <- all_pits[from, ]
    for (method in colnames(exact)) checks[[length(checks) + 1L]] <- data.frame(
      check = paste0(scenarios$scenario[matching[j]], "_analytic_", method),
      error = max(abs(all_pits[from, method] - exact[from - n_reference, method])),
      tolerance = reference_tolerance)
  }
}
quadrature <- pracma::gaussLegendre(64L, 0, 1)
t_nodes <- (quadrature$x / (sqrt(d) * (1 - quadrature$x)))^2
central_energy <- sqrt(d / pi) * sum(
  -expm1(-d/2 * log1p(2 * t_nodes)) * quadrature$w / quadrature$x^2)
checks[[length(checks) + 1L]] <- data.frame(check = "energy_exact_central_chi",
  error = abs(central_energy - sqrt(2) * exp(lgamma((d + 1)/2) - lgamma(d/2))),
  tolerance = 1e-10)
checks <- do.call(rbind, checks)
checks$passed <- is.finite(checks$error) & checks$error <= checks$tolerance
stopifnot(all(checks$passed))
pits <- data.frame(scenario = rep(scenarios$scenario, each = n),
                   case = rep(seq_len(n), nrow(scenarios)), pit_values)
stopifnot(all(is.finite(pit_values)), all(pit_values >= 0), all(pit_values <= 1))
settings <- data.frame(n_observations = n, n_reference = n_reference, n_null = n_null,
  seed = 20261015L, dimension = d, cdf_nodes = 48L, energy_nodes = 64L, bins = 20L,
  histogram_scale = "density", reference_cdf_tolerance = reference_tolerance,
  null_cdf_tolerance = null_tolerance)
write.csv(pits, file.path(result_dir, "population_pit_values.csv"), row.names = FALSE)
write.csv(checks, file.path(result_dir, "population_pit_checks.csv"), row.names = FALSE)
write.csv(scenarios, file.path(result_dir, "population_pit_scenarios.csv"), row.names = FALSE)
write.csv(settings, file.path(result_dir, "population_pit_settings.csv"), row.names = FALSE)
write.csv(do.call(rbind, null_tables), file.path(result_dir, "population_pit_null_values.csv"),
          row.names = FALSE)

# Pair low/high axis limits for each method and misspecification.
breaks <- seq(0, 1, length.out = 21L)
panels <- list()
summaries <- list()
pair_limits <- numeric()
for (scenario in scenarios$scenario) for (method in methods) {
  h <- hist(pits[pits$scenario == scenario, method], breaks = breaks,
            plot = FALSE, include.lowest = TRUE, right = TRUE)
  panels[[paste(scenario, method)]] <- h
  summaries[[length(summaries) + 1L]] <- data.frame(scenario = scenario, method = method,
    bin_left = head(breaks, -1), bin_right = tail(breaks, -1), count = h$counts,
    relative_frequency = h$counts / sum(h$counts), density = h$density)
  key <- paste(sub("_(low|high)$", "", scenario), method)
  previous <- if (key %in% names(pair_limits)) pair_limits[[key]] else 1.5
  pair_limits[key] <- max(previous, ceiling(max(h$density) / .5) * .5)
}
summaries <- do.call(rbind, summaries)
pair_key <- paste(sub("_(low|high)$", "", summaries$scenario), summaries$method)
summaries$paired_density_max <- unname(pair_limits[pair_key])
summaries$paired_relative_frequency_max <- summaries$paired_density_max / 20
summaries$uniform_reference_frequency <- 1/20
summaries$uniform_reference_density <- 1
write.csv(summaries, file.path(result_dir, "population_pit_histograms.csv"), row.names = FALSE)
for (scenario in scenarios$scenario) for (method in methods) {
  h <- panels[[paste(scenario, method)]]
  ymax <- pair_limits[[paste(sub("_(low|high)$", "", scenario), method)]]
  pdf(file.path(figure_dir, paste0("population_pit_", method, "_", scenario, ".pdf")),
      width = 72/72, height = 68/72, paper = "special", useDingbats = FALSE,
      title = paste("Pre-rank PIT:", method, scenario), compress = TRUE)
  par(mai = rep(0, 4), xaxs = "i", yaxs = "i", xpd = NA, lend = "butt", ps = 8)
  plot.new()
  plot.window(xlim = c(0, 72), ylim = c(0, 68))
  left <- 17; bottom <- 12; width <- 53; height <- 52
  rect(left + head(h$breaks, -1) * width, bottom,
       left + tail(h$breaks, -1) * width, bottom + h$density / ymax * height,
       col = "#167D8D", border = "white", lwd = .2)
  segments(left, bottom + height / ymax, left + width, bottom + height / ymax,
           col = "#C96539", lwd = .8, lty = "44")
  rect(left, bottom, left + width, bottom + height, border = "#19324A", col = NA, lwd = .5)
  xticks <- c(0, .5, 1)
  segments(left + xticks * width, bottom, left + xticks * width, bottom - 2,
           col = "#19324A", lwd = .5)
  text(left + xticks * width, bottom - 6, c("0", ".5", "1"), col = "#19324A", cex = .9)
  yticks <- unique(c(0, 1, ymax))
  segments(left - 2, bottom + yticks / ymax * height, left, bottom + yticks / ymax * height,
           col = "#19324A", lwd = .5)
  tick_labels <- vapply(yticks, format, character(1), trim = TRUE, scientific = FALSE)
  text(left - 3, bottom + yticks / ymax * height, tick_labels,
       adj = c(1, .5), col = "#19324A", cex = .85)
  dev.off()
}

# 6. Rosenblatt conditional PITs ------------------------------------------

rho_forecast <- .5
sigma_forecast <- matrix(c(1, rho_forecast, rho_forecast, 1), 2L)
sigma_true <- matrix(c(1, 1.25, 1.25, 1.75), 2L)
conditional_sd <- sqrt(1 - rho_forecast^2)
n <- 10000L
set.seed(20261008L)
observed <- mvtnorm::rmvnorm(n, sigma = sigma_true)
ink <- "#19324A"
observed_colors <- c("#167D8D", "#65AEB7")
forecast_colors <- c("#C96539", "#E6A069")
pale <- "#E3E8EC"

# Columns stay in original-coordinate order. Reversing the transform order
# reveals the misspecification hidden by the two uniform forward marginals.
observed_pits <- list(
  cbind(pnorm(observed[, 1L]),
        pnorm((observed[, 2L] - rho_forecast * observed[, 1L]) / conditional_sd)),
  cbind(pnorm((observed[, 1L] - rho_forecast * observed[, 2L]) / conditional_sd),
        pnorm(observed[, 2L])))
score_maps <- list(
  matrix(c(1, 0, -rho_forecast / conditional_sd, 1 / conditional_sd), 2L, byrow = TRUE),
  matrix(c(1 / conditional_sd, -rho_forecast / conditional_sd, 0, 1), 2L, byrow = TRUE))
breaks <- seq(0, 1, length.out = 11L)
histograms <- list()
score_covariances <- list()
histogram_values <- list()
maximum_density <- 0
bin_mass_error <- 0
bin_symmetry_error <- 0
reference_covariance_error <- 0
for (order in 1:2) {
  a <- score_maps[[order]]
  score_covariances[[order]] <- a %*% sigma_true %*% t(a)
  reference_covariance_error <- max(reference_covariance_error,
    abs(a %*% sigma_forecast %*% t(a) - diag(2L)))
  histograms[[order]] <- list()
  for (j in 1:2) {
    h <- hist(observed_pits[[order]][, j], breaks = breaks,
              plot = FALSE, include.lowest = TRUE)
    expected <- diff(pnorm(qnorm(breaks) / sqrt(score_covariances[[order]][j, j])))
    stopifnot(sum(h$counts) == n)
    bin_mass_error <- max(bin_mass_error, abs(sum(expected) - 1))
    bin_symmetry_error <- max(bin_symmetry_error, abs(expected - rev(expected)))
    histograms[[order]][[j]] <- h
    maximum_density <- max(maximum_density, h$density)
    step <- if (order == 1L) j else 3L - j
    histogram_values[[length(histogram_values) + 1L]] <- data.frame(
      order = c("1,2", "2,1")[order], coordinate = j, step = step,
      type = if (step == 2L) "conditional" else "marginal", sample_law = "Q", transform_law = "F",
      rho_true = cov2cor(sigma_true)[1L, 2L], rho_forecast = rho_forecast,
      true_variance_1 = sigma_true[1L, 1L], true_variance_2 = sigma_true[2L, 2L],
      true_covariance = sigma_true[1L, 2L],
      normal_score_variance = score_covariances[[order]][j, j],
      n_observed = n, seed = 20261008L, bin_left = head(breaks, -1L),
      bin_right = tail(breaks, -1L), count = h$counts, empirical_probability = h$counts / n,
      density = h$density, exact_probability = expected, exact_density = expected / diff(breaks))
  }
}
histogram_values <- do.call(rbind, histogram_values)
checks <- data.frame(
  check = c("forward_score_variances_equal_one", "reverse_score_variances_identity_error",
    "forward_score_correlation_identity_error", "pit_bin_mass_error", "pit_bin_symmetry_error",
    "forecast_reference_identity_covariance_error", "observed_bins_max_probability_error",
    "forward_pit_correlation_mc_error"),
  value = c(max(abs(diag(score_covariances[[1L]]) - 1)),
    max(abs(diag(score_covariances[[2L]]) - c(.25, 1.75))),
    abs(cov2cor(score_covariances[[1L]])[1L, 2L] - sqrt(3)/2), bin_mass_error,
    bin_symmetry_error, reference_covariance_error,
    max(abs(histogram_values$empirical_probability - histogram_values$exact_probability)),
    abs(cor(observed_pits[[1L]])[1L, 2L] - 6/pi * asin(sqrt(3)/4))),
  tolerance = c(rep(1e-12, 6L), .025, .025))
checks$passed <- checks$value <= checks$tolerance
stopifnot(all(checks$passed))
write.csv(histogram_values, file.path(result_dir, "rosenblatt_order_histograms.csv"), row.names = FALSE)
write.csv(checks, file.path(result_dir, "rosenblatt_order_checks.csv"), row.names = FALSE)

grDevices::cairo_pdf(file.path(figure_dir, "rosenblatt_order.pdf"),
                    width = 13.8, height = 2.4, pointsize = 16,
                    family = "sans", bg = "white")
layout(matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7),
              nrow = 3L, byrow = TRUE), heights = c(.13, .13, .74))
par(family = "sans", cex = 1, col = ink, col.axis = ink, col.lab = ink,
    col.main = ink, bg = "white", las = 1, lend = "round",
    mgp = c(1.15, .25, 0), tcl = -.18, xaxs = "i", yaxs = "i")
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", horiz = TRUE, bty = "n", lwd = c(7, 2.5), lty = c(1, 2),
       col = c(observed_colors[1], forecast_colors[1]),
       seg.len = 1.9, x.intersp = .65, text.col = ink, cex = .83,
       legend = c("Observed PITs (n = 10,000)", "Forecast reference: uniform"))
for (order in 1:2) {
  plot.new()
  text(.5, .45, if (order == 1L) expression("Order " * 1 %->% 2) else
         expression("Order " * 2 %->% 1), font = 2, cex = .95)
}
pit_labels <- expression(F[1](Y[1]), F[2*"|"*1](Y[2]*"|"*Y[1]),
                         F[2](Y[2]), F[1*"|"*2](Y[1]*"|"*Y[2]))
ymax <- max(2.5, ceiling(maximum_density * 12) / 10)
for (order in 1:2) for (step in 1:2) {
  j <- if (order == 1L) step else 3L - step
  panel <- 2L * (order - 1L) + step
  h <- histograms[[order]][[j]]
  par(mai = c(.61, .60, .025, .12))
  plot(NA, xlim = c(0, 1), ylim = c(0, ymax), axes = FALSE,
       xlab = "", ylab = "Density", cex.lab = .82)
  abline(h = c(0, 1, 2), col = pale, lwd = .7)
  rect(head(breaks, -1L), 0, tail(breaks, -1L), h$density,
       col = adjustcolor(observed_colors[order], alpha.f = .30),
       border = observed_colors[order], lwd = .85)
  abline(h = 1, col = forecast_colors[order], lty = 2, lwd = 2.2)
  text(.5, .90 * ymax, if (step == 1L) "Step 1: marginal" else
         "Step 2: conditional", cex = .76, col = ink)
  axis(1, at = c(0, .5, 1), labels = c("0", ".5", "1"), cex.axis = .77)
  axis(2, at = 0:2, cex.axis = .77)
  box(bty = "l", col = ink, lwd = .8)
  mtext(pit_labels[panel], side = 1, line = 1.45, cex = .87)
}
dev.off()

message("Done: slide PDFs in figures/; numerical checks in output/results/.")

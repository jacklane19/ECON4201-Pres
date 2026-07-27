# =============================================================================
# Three-variable proxy SVAR for Bishop & Tulip (2017)
# =============================================================================
#
# Reduced-form VAR:
#
#   x_t = c + A_1 x_(t-1) + ... + A_p x_(t-p) + u_t,
#
#   x_t = [log(real GDP), quarterly underlying inflation, cash-rate level]'.
#
# External instrument:
#
#   z_t = the raw residual from Equation (9), supplied in task1_shocks.csv.
#
# The proxy is NOT included as a variable in the VAR. Instead, it identifies the
# monetary-policy column of the structural impact matrix from
#
#   E[u_t z_t] = b_MP E[epsilon_MP,t z_t].
#
# The impact vector is normalised so that the cash rate rises by 1 percentage
# point on impact. GDP and inflation are unrestricted on impact in the baseline.
#
# Baseline choices:
#   - 3 variables: log real GDP, quarterly underlying inflation, cash rate;
#   - cash-rate measure: cr_qtr, the quarterly cash-rate series used in the
#     paper's conventional VAR;
#   - 3 VAR lags and a constant;
#   - sample: 1991Q2-2015Q4;
#   - raw Equation (9) residual as the external instrument;
#   - 20-quarter IRFs;
#   - 90% pointwise percentile intervals from a paired residual/proxy bootstrap.
#
# Bootstrap interpretation:
#   Each draw resamples the estimated VAR residual vector jointly with the proxy,
#   simulates a new VAR sample, re-estimates the VAR, and re-identifies the impact
#   vector. This captures VAR and proxy-moment uncertainty, but treats the supplied
#   Equation (9) residual series as the observed proxy. It does not re-estimate
#   Equation (9) inside every bootstrap draw.
# =============================================================================
library(readxl)
library(dplyr)
library(tidyr)
library(vars)
library(sandwich)
library(lmtest)
library(ggplot2)
library(patchwork)
# Namespace note: dplyr/tidyr verbs are qualified below because packages such as
# MASS can mask dplyr::select() in some R sessions.
# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------
data_file  <- "data.xlsx"

# All generated files are written to an "output" subfolder; task1_shocks.csv is
# read from there, having been produced by Shock_creation.R.
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outp <- function(f) file.path(out_dir, f)
shock_file <- outp("task1_shocks.csv")

p_lags   <- 3
H        <- 20
ci_level <- 0.90
boot_reps <- 1000       # use 199 while testing the script, then return to 1000
seed      <- 123
sample_start <- 1991.25  # 1991Q2
sample_end   <- 2015.75  # 2015Q4
# Quarterly cash-rate level. Change this to "cr_new" for a timing robustness.
policy_var <- "cr_qtr"
# Reject explosive bootstrap VARs. This conditions inference on a stable VAR.
enforce_stability <- TRUE
max_boot_attempts <- boot_reps * 8
stopifnot(file.exists(data_file), file.exists(shock_file))
set.seed(seed)
# -----------------------------------------------------------------------------
# 2. Import and transform the data
# -----------------------------------------------------------------------------
q <- as.data.frame(read_excel(data_file, sheet = "STATA_q"))
shocks <- read.csv(shock_file)
required_q <- c("year", "quarter", "rgdp", "cpi_trim", policy_var)
required_z <- c("year", "quarter", "shock")
if (!all(required_q %in% names(q))) {
  stop("The STATA_q sheet is missing one or more required variables: ",
       paste(setdiff(required_q, names(q)), collapse = ", "))
}
if (!all(required_z %in% names(shocks))) {
  stop("The shock file is missing one or more required variables: ",
       paste(setdiff(required_z, names(shocks)), collapse = ", "))
}
if (anyDuplicated(shocks[, c("year", "quarter")])) {
  stop("The shock file contains duplicate year-quarter observations.")
}
shocks <- shocks[, required_z]
names(shocks)[names(shocks) == "shock"] <- "proxy"
# Construct inflation before restricting the sample so that 1991Q2 inflation
# correctly uses the 1991Q1 price level.
d <- q %>%
  dplyr::arrange(year, quarter) %>%
  dplyr::mutate(
    qidx    = year + (quarter - 1) / 4,
    ln_rgdp = log(rgdp),
    infl    = log(cpi_trim) - dplyr::lag(log(cpi_trim))
  ) %>%
  dplyr::left_join(shocks, by = c("year", "quarter"))
var_names <- c("ln_rgdp", "infl", policy_var)
var_df <- d %>%
  dplyr::filter(qidx >= sample_start, qidx <= sample_end) %>%
  dplyr::filter(complete.cases(dplyr::across(dplyr::all_of(c(var_names, "proxy"))))) %>%
  dplyr::select(year, quarter, qidx, dplyr::all_of(var_names), proxy)
if (nrow(var_df) != 99) {
  warning("Expected 99 quarterly observations from 1991Q2 to 2015Q4, but found ",
          nrow(var_df), ". Check the data and sample alignment.")
}
if (nrow(var_df) <= p_lags + 10) {
  stop("Too few observations to estimate the requested VAR.")
}
cat(sprintf(
  "VAR sample: %dQ%d-%dQ%d; T = %d; p = %d\n",
  var_df$year[1], var_df$quarter[1],
  var_df$year[nrow(var_df)], var_df$quarter[nrow(var_df)],
  nrow(var_df), p_lags
))
# -----------------------------------------------------------------------------
# 3. Estimate the reduced-form VAR
# -----------------------------------------------------------------------------
Y_full <- as.matrix(var_df[, var_names])
colnames(Y_full) <- var_names
var_ts <- ts(
  Y_full,
  start = c(var_df$year[1], var_df$quarter[1]),
  frequency = 4
)
fit_var <- VAR(var_ts, p = p_lags, type = "const")
# Reduced-form innovations u_t. The VAR loses the first p observations.
U <- residuals(fit_var)
U <- as.matrix(U)
colnames(U) <- var_names
proxy_used <- var_df$proxy[(p_lags + 1):nrow(var_df)]
proxy_dates <- var_df[(p_lags + 1):nrow(var_df), c("year", "quarter", "qidx")]
stopifnot(nrow(U) == length(proxy_used))
stopifnot(all(complete.cases(U)), all(is.finite(proxy_used)))
root_modulus <- roots(fit_var, modulus = TRUE)
max_root <- max(root_modulus)
cat(sprintf("Largest companion-matrix root modulus: %.4f\n", max_root))
if (max_root >= 1) {
  warning("The estimated VAR is not stable. Long-horizon IRFs may be explosive.")
}
# -----------------------------------------------------------------------------
# 4. Proxy relevance: first stage using the cash-rate innovation
# -----------------------------------------------------------------------------
first_stage_data <- data.frame(
  year = proxy_dates$year,
  quarter = proxy_dates$quarter,
  proxy = proxy_used,
  policy_innovation = U[, policy_var]
)
first_stage <- lm(policy_innovation ~ proxy, data = first_stage_data)
# The VAR residual should be close to serially uncorrelated, but a one-lag
# Newey-West covariance is reported as a conservative diagnostic.
first_stage_vcov <- NeweyWest(
  first_stage,
  lag = 1,
  prewhite = FALSE,
  adjust = TRUE
)
first_stage_test <- coeftest(first_stage, vcov. = first_stage_vcov)
proxy_slope <- unname(first_stage_test["proxy", "Estimate"])
proxy_se    <- unname(first_stage_test["proxy", "Std. Error"])
proxy_t     <- unname(first_stage_test["proxy", "t value"])
robust_F    <- proxy_t^2
partial_R2  <- summary(first_stage)$r.squared
first_stage_summary <- data.frame(
  observations = nobs(first_stage),
  slope = proxy_slope,
  hac_standard_error = proxy_se,
  robust_t = proxy_t,
  robust_F = robust_F,
  partial_R_squared = partial_R2,
  correlation = cor(first_stage_data$proxy, first_stage_data$policy_innovation)
)
cat("\n=== Proxy first-stage diagnostic ===\n")
print(first_stage_test)
cat(sprintf("Robust first-stage F = %.2f\n", robust_F))
cat(sprintf("Partial R-squared     = %.3f\n", partial_R2))
if (robust_F < 10) {
  warning("The proxy first stage is weak by the conventional F=10 benchmark. ",
          "Treat normalized IRFs and ordinary bootstrap intervals cautiously.")
}
write.csv(first_stage_summary, outp("proxy_svar_first_stage.csv"), row.names = FALSE)
write.csv(first_stage_data, outp("proxy_svar_first_stage_data.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------
# 5. Identify the monetary-policy impact vector
# -----------------------------------------------------------------------------
proxy_moment <- function(residual_matrix, instrument) {
  # Sample analogue of E[u_t z_t], after demeaning both objects.
  U_c <- sweep(residual_matrix, 2, colMeans(residual_matrix), FUN = "-")
  z_c <- instrument - mean(instrument)
  drop(crossprod(U_c, z_c) / (length(z_c) - 1))
}
moment_hat <- proxy_moment(U, proxy_used)
names(moment_hat) <- var_names
if (abs(moment_hat[policy_var]) < 1e-10) {
  stop("The covariance between the proxy and cash-rate innovation is too close ",
       "to zero to normalize the proxy SVAR.")
}
# Only the relative scale of the impact vector is identified by the proxy.
# Normalize the cash-rate response to +1 percentage point on impact.
impact_vector <- moment_hat / moment_hat[policy_var]
names(impact_vector) <- var_names
proxy_correlations <- sapply(var_names, function(v) cor(U[, v], proxy_used))
impact_results <- data.frame(
  variable = var_names,
  proxy_covariance = unname(moment_hat[var_names]),
  proxy_correlation = unname(proxy_correlations[var_names]),
  normalized_impact = unname(impact_vector[var_names])
)
cat("\n=== Proxy moments and normalized impact vector ===\n")
print(impact_results, digits = 5, row.names = FALSE)
write.csv(impact_results, outp("proxy_svar_impact_vector.csv"), row.names = FALSE)
# -----------------------------------------------------------------------------
# 6. Propagate the identified impact vector through the VAR
# -----------------------------------------------------------------------------
make_structural_irfs <- function(var_fit, impact, horizons) {
  Phi_array <- Phi(var_fit, nstep = horizons)
  K <- length(impact)
  
  raw_irf <- matrix(
    NA_real_,
    nrow = horizons + 1,
    ncol = K,
    dimnames = list(horizon = 0:horizons, response = names(impact))
  )
  
  for (h in 0:horizons) {
    raw_irf[h + 1, ] <- Phi_array[, , h + 1] %*% impact
  }
  
  raw_irf
}
format_irfs <- function(raw_irf, policy_name) {
  # ln_rgdp is a log level, so multiply its response by 100 to report percent.
  # infl is a raw quarterly log difference, so multiply by 100 to report ppt.
  inflation_ppt <- 100 * raw_irf[, "infl"]
  
  data.frame(
    horizon = 0:(nrow(raw_irf) - 1),
    real_gdp = 100 * raw_irf[, "ln_rgdp"],
    inflation = inflation_ppt,
    price_level = cumsum(inflation_ppt),
    cash_rate = raw_irf[, policy_name]
  )
}
raw_irf_hat <- make_structural_irfs(fit_var, impact_vector, H)
point_wide <- format_irfs(raw_irf_hat, policy_var)
response_names <- c("real_gdp", "inflation", "price_level", "cash_rate")
response_labels <- c(
  real_gdp = "Real GDP (% deviation)",
  inflation = "Underlying inflation (quarterly ppt)",
  price_level = "Price level (% deviation)",
  cash_rate = "Cash rate (ppt)"
)
# -----------------------------------------------------------------------------
# 7. Paired residual/proxy bootstrap
# -----------------------------------------------------------------------------
# Extract estimated VAR dynamics for simulation.
A_hat <- Acoef(fit_var)
const_hat <- sapply(var_names, function(v) {
  unname(coef(fit_var$varresult[[v]])["const"])
})
names(const_hat) <- var_names
# Center the empirical innovations before resampling.
U_centered <- sweep(U, 2, colMeans(U), FUN = "-")
proxy_centered <- proxy_used - mean(proxy_used)
simulate_var_sample <- function(initial_values, A_list, intercept, innovations) {
  p <- length(A_list)
  T_total <- p + nrow(innovations)
  K <- ncol(initial_values)
  
  Y_sim <- matrix(NA_real_, nrow = T_total, ncol = K)
  colnames(Y_sim) <- colnames(initial_values)
  Y_sim[1:p, ] <- initial_values[1:p, , drop = FALSE]
  
  for (tt in (p + 1):T_total) {
    prediction <- as.numeric(intercept)
    for (lag_number in seq_len(p)) {
      prediction <- prediction + drop(
        A_list[[lag_number]] %*% Y_sim[tt - lag_number, ]
      )
    }
    Y_sim[tt, ] <- prediction + innovations[tt - p, ]
  }
  
  Y_sim
}
boot_irfs <- array(
  NA_real_,
  dim = c(boot_reps, H + 1, length(response_names)),
  dimnames = list(
    draw = seq_len(boot_reps),
    horizon = 0:H,
    response = response_names
  )
)
accepted <- 0L
attempts <- 0L
failed_unstable <- 0L
failed_weak <- 0L
failed_other <- 0L
N_resid <- nrow(U_centered)
initial_values <- Y_full[1:p_lags, , drop = FALSE]
cat(sprintf("\nStarting %d-draw paired residual/proxy bootstrap...\n", boot_reps))
while (accepted < boot_reps && attempts < max_boot_attempts) {
  attempts <- attempts + 1L
  
  # Resample the residual vector and proxy jointly. This preserves the empirical
  # relationship between the external instrument and all VAR innovations.
  draw_index <- sample.int(N_resid, size = N_resid, replace = TRUE)
  U_star <- U_centered[draw_index, , drop = FALSE]
  z_star <- proxy_centered[draw_index]
  
  Y_star <- simulate_var_sample(
    initial_values = initial_values,
    A_list = A_hat,
    intercept = const_hat,
    innovations = U_star
  )
  
  fit_star <- try(
    VAR(
      ts(Y_star, start = c(var_df$year[1], var_df$quarter[1]), frequency = 4),
      p = p_lags,
      type = "const"
    ),
    silent = TRUE
  )
  
  if (inherits(fit_star, "try-error")) {
    failed_other <- failed_other + 1L
    next
  }
  
  if (enforce_stability && max(roots(fit_star, modulus = TRUE)) >= 1) {
    failed_unstable <- failed_unstable + 1L
    next
  }
  
  U_hat_star <- as.matrix(residuals(fit_star))
  colnames(U_hat_star) <- var_names
  
  if (nrow(U_hat_star) != length(z_star)) {
    failed_other <- failed_other + 1L
    next
  }
  
  moment_star <- proxy_moment(U_hat_star, z_star)
  names(moment_star) <- var_names
  
  # Skip only genuinely unidentified/non-finite draws. Small denominators are
  # otherwise retained, because extreme IRFs are informative about weak-proxy
  # uncertainty rather than a numerical problem.
  if (!all(is.finite(moment_star)) ||
      abs(moment_star[policy_var]) < .Machine$double.eps^0.5) {
    failed_weak <- failed_weak + 1L
    next
  }
  
  impact_star <- moment_star / moment_star[policy_var]
  names(impact_star) <- var_names
  
  raw_star <- try(make_structural_irfs(fit_star, impact_star, H), silent = TRUE)
  if (inherits(raw_star, "try-error") || !all(is.finite(raw_star))) {
    failed_other <- failed_other + 1L
    next
  }
  
  formatted_star <- format_irfs(raw_star, policy_var)
  
  accepted <- accepted + 1L
  boot_irfs[accepted, , ] <- as.matrix(formatted_star[, response_names])
  
  if (accepted %% 100 == 0 || accepted == boot_reps) {
    cat(sprintf("  accepted %d of %d draws\n", accepted, boot_reps))
  }
}
if (accepted < boot_reps) {
  warning("Only ", accepted, " bootstrap draws were accepted after ", attempts,
          " attempts. Confidence intervals use the accepted draws.")
}
if (accepted < 100) {
  stop("Too few successful bootstrap draws for meaningful confidence intervals.")
}
boot_irfs <- boot_irfs[seq_len(accepted), , , drop = FALSE]
alpha <- 1 - ci_level
band_list <- lapply(response_names, function(response_name) {
  response_draws <- boot_irfs[, , response_name, drop = FALSE]
  response_draws <- matrix(response_draws, nrow = accepted, ncol = H + 1)
  
  data.frame(
    horizon = 0:H,
    response = response_name,
    lower = apply(
      response_draws,
      2,
      quantile,
      probs = alpha / 2,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ),
    upper = apply(
      response_draws,
      2,
      quantile,
      probs = 1 - alpha / 2,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
})
bands_long <- dplyr::bind_rows(band_list)
point_long <- point_wide %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(response_names),
    names_to = "response",
    values_to = "estimate"
  )
irf_results <- point_long %>%
  dplyr::left_join(bands_long, by = c("horizon", "response")) %>%
  dplyr::mutate(
    response_label = unname(response_labels[response]),
    response_label = factor(
      response_label,
      levels = unname(response_labels[response_names])
    )
  ) %>%
  dplyr::arrange(factor(response, levels = response_names), horizon)
bootstrap_summary <- data.frame(
  requested_draws = boot_reps,
  accepted_draws = accepted,
  attempts = attempts,
  rejected_unstable = failed_unstable,
  rejected_near_zero_proxy_moment = failed_weak,
  rejected_other = failed_other,
  stability_enforced = enforce_stability,
  seed = seed
)
write.csv(irf_results, outp("proxy_svar_irfs.csv"), row.names = FALSE)
write.csv(bootstrap_summary, outp("proxy_svar_bootstrap_summary.csv"), row.names = FALSE)
key_horizons <- irf_results %>%
  dplyr::filter(horizon %in% c(4, 8)) %>%
  dplyr::select(response, response_label, horizon, estimate, lower, upper)
write.csv(key_horizons, outp("proxy_svar_h4_h8.csv"), row.names = FALSE)
cat("\n=== Proxy-SVAR responses at horizons 4 and 8 ===\n")
print(key_horizons, digits = 4, row.names = FALSE)
# -----------------------------------------------------------------------------
# 8. Plots
# -----------------------------------------------------------------------------
# Per-response accent palette (keys must match response_labels exactly).
irf_palette <- c(
  "Real GDP (% deviation)"               = "#2C6E9B",
  "Underlying inflation (quarterly ppt)" = "#C1553B",
  "Price level (% deviation)"            = "#B23A48",
  "Cash rate (ppt)"                      = "#4E6E58"
)
# Shared clean theme.
theme_svar <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, colour = "#16181c"),
    plot.subtitle    = element_text(size = 9.5, colour = "#5f6368", margin = margin(b = 6)),
    plot.caption     = element_text(size = 8, colour = "#9aa0a6", hjust = 0),
    strip.text       = element_text(face = "bold", size = 11, colour = "#20232a", hjust = 0),
    axis.title       = element_text(size = 8.5, colour = "#5f6368"),
    axis.text        = element_text(size = 8.5, colour = "#5f6368"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#eef0f2", linewidth = 0.5),
    axis.ticks       = element_blank(),
    legend.position  = "none",
    panel.spacing    = unit(1.1, "lines"),
    plot.margin      = margin(8, 14, 8, 8)
  )
# Panel display order: cash rate, GDP, inflation, price level.
irf_results$response_label <- factor(
  irf_results$response_label,
  levels = c("Cash rate (ppt)", "Real GDP (% deviation)",
             "Underlying inflation (quarterly ppt)", "Price level (% deviation)")
)
irf_markers <- subset(irf_results, horizon %in% c(4, 8))
irf_plot <- ggplot(irf_results,
                   aes(x = horizon, y = estimate,
                       colour = response_label, fill = response_label)) +
  geom_hline(yintercept = 0, colour = "#9aa0a6", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(data = irf_markers, shape = 21, colour = "white", size = 2.3, stroke = 0.7) +
  facet_wrap(~ response_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = irf_palette, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, H, by = 4), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Proxy-SVAR responses to a contractionary monetary policy shock",
    subtitle = paste0(
      "External instrument: Equation (9) residual; VAR = [log GDP, inflation, ",
      policy_var, "]. Impact normalised to +1 ppt cash rate; shaded = ",
      round(100 * ci_level), "% paired-bootstrap band."
    ),
    x = "Quarters after shock",
    y = NULL,
    caption = paste0(
      "Sample 1991Q2-2015Q4; VAR lags: ", p_lags,
      ". Price level = cumulative inflation response. Dots mark h = 4 and 8."
    )
  ) +
  theme_svar
first_stage_plot <- ggplot(first_stage_data, aes(x = proxy, y = policy_innovation)) +
  geom_hline(yintercept = 0, colour = "#c9ccd1", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#c9ccd1", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#B23A48", linewidth = 1.1) +
  geom_point(shape = 21, fill = "#4E6E58", colour = "white",
             size = 2.3, stroke = 0.4, alpha = 0.85) +
  labs(
    title = "Proxy relevance",
    subtitle = sprintf(
      "Cash-rate VAR innovation on the Equation (9) residual   |   robust F = %.1f, partial R² = %.2f",
      robust_F, partial_R2
    ),
    x = "Equation (9) residual (ppt)",
    y = paste0(policy_var, " reduced-form innovation (ppt)")
  ) +
  theme_svar
combined_plot <- irf_plot / first_stage_plot +
  plot_layout(heights = c(2.2, 1))
ggsave(outp("extension_proxy_svar_irfs.png"), irf_plot,         width = 10.5, height = 7,  dpi = 180, bg = "white")
ggsave(outp("proxy_svar_first_stage.png"),    first_stage_plot, width = 7.5,  height = 5,  dpi = 180, bg = "white")
ggsave(outp("proxy_svar_results.png"),        combined_plot,    width = 10.5, height = 10, dpi = 180, bg = "white")
# -----------------------------------------------------------------------------
# 9. Save all objects for later comparison with the direct-shock VAR and LP
# -----------------------------------------------------------------------------
saveRDS(
  list(
    settings = list(
      p_lags = p_lags,
      H = H,
      ci_level = ci_level,
      boot_reps_requested = boot_reps,
      policy_var = policy_var,
      sample_start = sample_start,
      sample_end = sample_end,
      seed = seed
    ),
    data = var_df,
    fit_var = fit_var,
    reduced_form_residuals = U,
    proxy_used = proxy_used,
    first_stage = first_stage,
    first_stage_vcov = first_stage_vcov,
    first_stage_summary = first_stage_summary,
    proxy_moment = moment_hat,
    impact_vector = impact_vector,
    raw_irf = raw_irf_hat,
    irf_results = irf_results,
    bootstrap_irfs = boot_irfs,
    bootstrap_summary = bootstrap_summary,
    root_modulus = root_modulus
  ),
  outp("proxy_svar_results.rds")
)
cat("\nSaved to", out_dir, ":\n")
cat("  proxy_svar_first_stage.csv\n")
cat("  proxy_svar_first_stage_data.csv\n")
cat("  proxy_svar_impact_vector.csv\n")
cat("  proxy_svar_irfs.csv\n")
cat("  proxy_svar_h4_h8.csv\n")
cat("  proxy_svar_bootstrap_summary.csv\n")
cat("  extension_proxy_svar_irfs.png\n")
cat("  proxy_svar_first_stage.png\n")
cat("  proxy_svar_results.png\n")
cat("  proxy_svar_results.rds\n")
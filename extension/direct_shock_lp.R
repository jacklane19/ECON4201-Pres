# Direct-shock local projections using the Bishop & Tulip (2017)
# Equation (9) monetary-policy shock series.
#
# Baseline specification:
#   - Treatment: raw Equation (9) residual, mpshock_t, measured in percentage points.
#   - Controls: 3 lags of GDP growth, quarterly inflation, cash-rate level,
#     and the monetary-policy shock.
#   - Inference: horizon-specific Newey-West standard errors with lag h + 1.
#   - Outcomes:
#       1. quarterly underlying inflation at t+h;
#       2. price-level change from t to t+h (h >= 1);
#       3. real-GDP level change from t-1 to t+h;
#       4. cash-rate change from t-1 to t+h.
#
# The price-level outcome follows the timing discussed in the accompanying notes:
# it excludes quarter-t inflation because the forecast-month policy decision occurs
# within quarter t. Its horizon-zero response is therefore normalised to zero.
#
# Important: the HAC intervals below treat the supplied shock series as observed.
# They do not incorporate uncertainty from estimating Equation (9). A block
# bootstrap that re-estimates Equation (9) is the natural robustness extension.

library(readxl)
library(dplyr)
library(sandwich)
library(lmtest)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------
data_file  <- "data.xlsx"

# All generated files are written to an "output" subfolder (task1_shocks.csv is
# read from there, having been produced by Shock_creation.R).
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outp <- function(f) file.path(out_dir, f)
shock_file <- outp("task1_shocks.csv")

p_lags   <- 3
H        <- 16
ci_level <- 0.90

shock_start <- 1991.25  # 1991Q2
shock_end   <- 2015.75  # 2015Q4

stopifnot(file.exists(data_file), file.exists(shock_file))

# -----------------------------------------------------------------------------
# 2. Read and transform the data
# -----------------------------------------------------------------------------
q <- as.data.frame(read_excel(data_file, sheet = "STATA_q"))
shocks <- read.csv(shock_file)

required_q <- c("year", "quarter", "rgdp", "cpi_trim", "cr_qtr")
required_s <- c("year", "quarter", "shock")
stopifnot(all(required_q %in% names(q)))
stopifnot(all(required_s %in% names(shocks)))
stopifnot(!anyDuplicated(shocks[, c("year", "quarter")]))

shocks <- shocks[, required_s]
names(shocks)[names(shocks) == "shock"] <- "mpshock"

d <- q %>%
  left_join(shocks, by = c("year", "quarter")) %>%
  arrange(year, quarter) %>%
  mutate(
    qidx       = year + (quarter - 1) / 4,
    ln_rgdp    = log(rgdp),
    ln_price   = log(cpi_trim),
    gdp_growth = 100 * (ln_rgdp - lag(ln_rgdp)),
    infl       = 100 * (ln_price - lag(ln_price)),
    in_shock_sample = qidx >= shock_start & qidx <= shock_end & !is.na(mpshock)
  )

# Construct predetermined controls W_{t-1}: three lags of each variable.
for (j in seq_len(p_lags)) {
  d[[paste0("gdp_growth_l", j)]] <- dplyr::lag(d$gdp_growth, j)
  d[[paste0("infl_l", j)]]       <- dplyr::lag(d$infl, j)
  d[[paste0("cr_qtr_l", j)]]     <- dplyr::lag(d$cr_qtr, j)
  d[[paste0("mpshock_l", j)]]    <- dplyr::lag(d$mpshock, j)
}

control_names <- unlist(lapply(seq_len(p_lags), function(j) {
  c(
    paste0("gdp_growth_l", j),
    paste0("infl_l", j),
    paste0("cr_qtr_l", j),
    paste0("mpshock_l", j)
  )
}))

# The first usable treatment quarter is normally 1992Q1 because three lagged
# shocks are included and the supplied shock series begins in 1991Q2.
base_complete <- d$in_shock_sample &
  complete.cases(d[, c("mpshock", control_names)])

cat(sprintf(
  "Potential treatment observations after lag controls: %d (%dQ%d to %dQ%d)\n",
  sum(base_complete),
  d$year[which(base_complete)[1]], d$quarter[which(base_complete)[1]],
  d$year[tail(which(base_complete), 1)], d$quarter[tail(which(base_complete), 1)]
))

# -----------------------------------------------------------------------------
# 3. Horizon-specific outcomes
# -----------------------------------------------------------------------------
make_lp_outcome <- function(data, horizon, outcome) {
  if (outcome == "inflation") {
    # Quarterly underlying inflation at t+h, in percentage points.
    return(dplyr::lead(data$infl, horizon))
  }
  if (outcome == "price_level") {
    # 100*[log(P_{t+h}) - log(P_t)] = sum of inflation from t+1 to t+h.
    # Horizon zero is a normalization and is handled outside the regression.
    if (horizon == 0) return(rep(0, nrow(data)))
    return(100 * (dplyr::lead(data$ln_price, horizon) - data$ln_price))
  }
  if (outcome == "real_gdp") {
    # GDP level at t+h relative to the quarter immediately before the shock.
    return(100 * (dplyr::lead(data$ln_rgdp, horizon) - dplyr::lag(data$ln_rgdp, 1)))
  }
  if (outcome == "cash_rate") {
    # Cash-rate path at t+h relative to the quarter immediately before the shock.
    return(dplyr::lead(data$cr_qtr, horizon) - dplyr::lag(data$cr_qtr, 1))
  }
  stop("Unknown outcome: ", outcome)
}

# Panel display order: cash rate, GDP, inflation, price level.
outcome_labels <- c(
  cash_rate   = "Cash rate (ppt from t-1)",
  real_gdp    = "Real GDP level (% from t-1)",
  inflation   = "Underlying inflation (quarterly ppt)",
  price_level = "Price level (% from quarter t)"
)

# -----------------------------------------------------------------------------
# 4. Estimate one direct-shock LP regression
# -----------------------------------------------------------------------------
estimate_one_lp <- function(data, horizon, outcome, ci = ci_level) {
  # The price-level response at h=0 is set to zero by timing convention rather
  # than estimated from a regression.
  if (outcome == "price_level" && horizon == 0) {
    return(data.frame(
      outcome = outcome,
      outcome_label = unname(outcome_labels[outcome]),
      horizon = horizon,
      estimate = 0,
      std_error = 0,
      lower = 0,
      upper = 0,
      n = NA_integer_,
      nw_lag = NA_integer_,
      r_squared = NA_real_,
      estimated = FALSE
    ))
  }
  
  work <- data
  work$lp_y <- make_lp_outcome(work, horizon, outcome)
  model_vars <- c("lp_y", "mpshock", control_names)
  keep <- work$in_shock_sample & complete.cases(work[, model_vars])
  est_data <- work[keep, model_vars, drop = FALSE]
  
  if (nrow(est_data) <= length(model_vars) + 2) {
    stop(sprintf(
      "Too few observations for outcome %s at horizon %d.", outcome, horizon
    ))
  }
  
  # LP_h: y_{t,h} = alpha_h + beta_h*mpshock_t + Gamma_h'W_{t-1} + u_{t,h}
  lp_formula <- reformulate(
    termlabels = c("mpshock", control_names),
    response = "lp_y"
  )
  fit <- lm(lp_formula, data = est_data)
  
  # Overlapping future outcomes induce serial correlation. h+1 is used as a
  # transparent horizon-specific HAC bandwidth; h is also a sensible robustness.
  nw_lag <- horizon + 1L
  vcov_nw <- NeweyWest(
    fit,
    lag = nw_lag,
    prewhite = FALSE,
    adjust = TRUE
  )
  test <- coeftest(fit, vcov. = vcov_nw)
  
  beta <- unname(test["mpshock", "Estimate"])
  se   <- unname(test["mpshock", "Std. Error"])
  critical <- qnorm(1 - (1 - ci) / 2)
  
  data.frame(
    outcome = outcome,
    outcome_label = unname(outcome_labels[outcome]),
    horizon = horizon,
    estimate = beta,
    std_error = se,
    lower = beta - critical * se,
    upper = beta + critical * se,
    n = nobs(fit),
    nw_lag = nw_lag,
    r_squared = summary(fit)$r.squared,
    estimated = TRUE
  )
}

# -----------------------------------------------------------------------------
# 5. Estimate the full impulse-response paths
# -----------------------------------------------------------------------------
outcomes <- names(outcome_labels)
results_list <- vector("list", length(outcomes) * (H + 1))
k <- 1L
for (outcome in outcomes) {
  for (h in 0:H) {
    results_list[[k]] <- estimate_one_lp(d, h, outcome)
    k <- k + 1L
  }
}

lp_results <- bind_rows(results_list) %>%
  arrange(factor(outcome, levels = outcomes), horizon)

write.csv(lp_results, outp("direct_shock_lp_results.csv"), row.names = FALSE)

key_horizons <- lp_results %>%
  filter(horizon %in% c(4, 8)) %>%
  select(outcome, outcome_label, horizon, estimate, std_error, lower, upper, n)
write.csv(key_horizons, outp("direct_shock_lp_h4_h8.csv"), row.names = FALSE)

cat("\n=== Direct-shock LP responses at horizons 4 and 8 ===\n")
print(key_horizons, digits = 4, row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Plot the responses
# -----------------------------------------------------------------------------
plot_data <- lp_results %>%
  mutate(
    outcome_label = factor(
      outcome_label,
      levels = unname(outcome_labels[outcomes])
    )
  )

# Per-outcome accent palette (keyed by the outcome_label factor levels).
lp_palette <- c(
  "Underlying inflation (quarterly ppt)" = "#C1553B",
  "Price level (% from quarter t)"       = "#B23A48",
  "Real GDP level (% from t-1)"          = "#2C6E9B",
  "Cash rate (ppt from t-1)"             = "#4E6E58"
)

# Shared clean theme.
theme_lp <- theme_minimal(base_size = 11) +
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

markers <- subset(plot_data, horizon %in% c(4, 8) & estimated)

p_all <- ggplot(plot_data, aes(horizon, estimate, colour = outcome_label, fill = outcome_label)) +
  geom_hline(yintercept = 0, colour = "#9aa0a6", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(data = markers, shape = 21, colour = "white", size = 2.3, stroke = 0.7) +
  facet_wrap(~ outcome_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = lp_palette, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, H, by = 4), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Direct-shock local projections",
    subtitle = paste0(
      "Response to a 1 percentage point Equation (9) monetary-policy shock; ",
      ci_level * 100, "% Newey-West intervals. Dots mark h = 4 and 8."
    ),
    x = "Quarters after shock",
    y = "Response",
    caption = paste0(
      "Controls: ", p_lags,
      " lags of GDP growth, inflation, cash rate and the policy shock. ",
      "HAC bandwidth = h + 1."
    )
  ) +
  theme_lp

ggsave(outp("extension_local_projections.png"), p_all, width = 10, height = 7, dpi = 160)

# Figure-9-comparable price-level path.
price_data <- filter(plot_data, outcome == "price_level")

p_price <- ggplot(price_data, aes(horizon, estimate)) +
  geom_hline(yintercept = 0, colour = "#9aa0a6", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#B23A48", alpha = 0.16) +
  geom_line(colour = "#B23A48", linewidth = 1.1, lineend = "round") +
  geom_point(data = filter(price_data, horizon %in% c(4, 8)),
             fill = "#B23A48", colour = "white", shape = 21, size = 2.7, stroke = 0.8) +
  scale_x_continuous(breaks = seq(0, H, by = 4), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Price-level response from direct-shock local projections",
    subtitle = paste0(
      "Cumulative inflation from t+1 to t+h after a 1 ppt policy shock; ",
      ci_level * 100, "% Newey-West intervals"
    ),
    x = "Quarters after shock",
    y = "% deviation from the quarter-t price level"
  ) +
  theme_lp

ggsave(outp("local_projections_price.png"), p_price, width = 8, height = 5.5, dpi = 160)

saveRDS(
  list(
    data = d,
    control_names = control_names,
    settings = list(
      p_lags = p_lags,
      H = H,
      ci_level = ci_level,
      shock_start = shock_start,
      shock_end = shock_end
    ),
    lp_results = lp_results,
    key_horizons = key_horizons,
    plots = list(all = p_all, price_level = p_price)
  ),
  outp("direct_shock_lp_results.rds")
)

cat("\nSaved to", out_dir, ":\n")
cat("  direct_shock_lp_results.csv\n")
cat("  direct_shock_lp_h4_h8.csv\n")
cat("  direct_shock_lp_all.png\n")
cat("  direct_shock_lp_price_level.png\n")
cat("  direct_shock_lp_results.rds\n")
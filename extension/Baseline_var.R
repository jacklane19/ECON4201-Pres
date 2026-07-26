# Bishop & Tulip (2017), Section 3.1 / Figure 9: baseline "b"
#
# IMPORTANT REPLICATION DETAIL
# ----------------------------
# The paper says that the policy-shock series is cumulated before entering the VAR.
# The public Stata file constructs *_c series, but its Figure 9 VAR loop actually
# estimates: var ln_rgdp infl mpshock, not mpshock_c. The RAW-shock specification
# below reproduces the published baseline point estimates (about 0.506 at h=4 and
# 0.571 at h=8). A separate paper-text specification is included at the end.

library(readxl)
library(vars)
library(ggplot2)

# All generated files are written to an "output" subfolder (task1_shocks.csv is
# read from there, having been produced by Shock_creation.R).
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outp <- function(f) file.path(out_dir, f)

set.seed(123)
q <- as.data.frame(read_excel("data.xlsx", sheet = "STATA_q"))
shocks <- read.csv(outp("task1_shocks.csv"))

q <- q[order(q$year, q$quarter), ]
q$ln_rgdp <- log(q$rgdp)
q$infl <- c(NA, diff(log(q$cpi_trim)))

merged <- merge(
  q, shocks[, c("year", "quarter", "shock")],
  by = c("year", "quarter"), all.x = TRUE
)
merged <- merged[order(merged$year, merged$quarter), ]
merged$qidx <- merged$year + (merged$quarter - 1) / 4
merged$mpshock <- merged$shock

# The released baseline shocks are gap-free from 1991Q2 to 2016Q1.
# Keep the zero-fill used by the broader Stata workflow, but assert that it is not
# doing substantive work for the supplied baseline series.
in_shock_window <- merged$qidx >= 1991.25 & merged$qidx <= 2016.00
stopifnot(!any(is.na(merged$mpshock[in_shock_window])))
merged$mpshock[in_shock_window & is.na(merged$mpshock)] <- 0

# Figure 9 sample: 1991Q2-2015Q4. With p=3, the first effective VAR observation
# is 1992Q1, as in Stata's L3.mpshock != . condition.
keep <- merged$qidx >= 1991.25 & merged$qidx <= 2015.75 &
  complete.cases(merged[, c("ln_rgdp", "infl", "mpshock")])
var_df <- merged[keep, c("year", "quarter", "ln_rgdp", "infl", "mpshock")]
stopifnot(nrow(var_df) == 99)
stopifnot(var_df$year[1] == 1991, var_df$quarter[1] == 2)
stopifnot(var_df$year[nrow(var_df)] == 2015, var_df$quarter[nrow(var_df)] == 4)

# Cholesky ordering: GDP, underlying inflation, policy shock.
var_ts <- ts(
  var_df[, c("ln_rgdp", "infl", "mpshock")],
  start = c(1991, 2), frequency = 4
)
fit_var <- VAR(var_ts, p = 3, type = "const")

# Stata uses step(30), bsp reps(1000). The vars package also provides a residual
# bootstrap, but its reported bands are percentile bands. Stata's do-file instead
# forms coirf +/- 1.645*stdcoirf. Thus point estimates below are directly comparable;
# the R confidence limits are close analogues, not bit-for-bit Stata replications.
h <- 30
impact_scale <- irf(
  fit_var, impulse = "mpshock", response = "mpshock",
  n.ahead = 1, ortho = TRUE, boot = FALSE
)$irf$mpshock[1, "mpshock"]

scaled_irf <- function(response, cumulative = FALSE, multiply_100 = TRUE,
                       n.ahead = h, runs = 1000, ci = 0.90, seed = 123) {
  x <- irf(
    fit_var, impulse = "mpshock", response = response,
    n.ahead = n.ahead, ortho = TRUE, cumulative = cumulative,
    boot = TRUE, runs = runs, ci = ci, seed = seed
  )
  multiplier <- (if (multiply_100) 100 else 1) / impact_scale
  x$irf$mpshock <- x$irf$mpshock * multiplier
  x$Lower$mpshock <- x$Lower$mpshock * multiplier
  x$Upper$mpshock <- x$Upper$mpshock * multiplier
  attr(x, "response") <- response
  x
}

# Figure 9 uses the PRICE LEVEL: cumulative orthogonalised response of inflation.
irf_price <- scaled_irf("infl", cumulative = TRUE, multiply_100 = TRUE)

# These additional objects distinguish economically meaningful level responses
# from the public do-file's generic use of coirf for every response variable.
irf_infl <- scaled_irf("infl", cumulative = FALSE, multiply_100 = TRUE)
irf_gdp_level <- scaled_irf("ln_rgdp", cumulative = FALSE, multiply_100 = TRUE)
irf_gdp_do_coirf <- scaled_irf("ln_rgdp", cumulative = TRUE, multiply_100 = TRUE)
irf_shock <- scaled_irf("mpshock", cumulative = FALSE, multiply_100 = FALSE)

tidy_irf <- function(x, response = attr(x, "response")) {
  data.frame(
    horizon = 0:(nrow(x$irf$mpshock) - 1),
    estimate = x$irf$mpshock[, response],
    lower = x$Lower$mpshock[, response],
    upper = x$Upper$mpshock[, response]
  )
}

price_path <- tidy_irf(irf_price, "infl")
figure9_b <- price_path[price_path$horizon %in% c(4, 8), ]
figure9_b$series <- "b: baseline"

cat("\n=== Figure 9 baseline b: price-level response ===\n")
print(figure9_b, digits = 6, row.names = FALSE)
cat("Expected point estimates from the released data: h4 about 0.506; h8 about 0.571.\n")
stopifnot(abs(figure9_b$estimate[figure9_b$horizon == 4] - 0.5057) < 0.01)
stopifnot(abs(figure9_b$estimate[figure9_b$horizon == 8] - 0.5713) < 0.01)
write.csv(figure9_b, outp("figure9_baseline_b.csv"), row.names = FALSE)

# Time path underlying the two baseline bars.
p_path <- ggplot(price_path, aes(horizon, estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18) +
  geom_line(linewidth = 0.8) +
  geom_point(data = figure9_b, size = 2) +
  labs(
    title = "Figure 9 baseline b: price-level response",
    subtitle = "Raw policy-shock series, VAR(3), Cholesky order GDP-inflation-shock",
    x = "Quarters after shock", y = "% deviation from baseline",
    caption = "Bands are vars-package 90% percentile bootstrap intervals; Stata uses +/-1.645 bootstrap SE."
  ) +
  theme_minimal(base_size = 11)
ggsave(outp("figure9_baseline_path_reviewed.png"), p_path,
       width = 8, height = 5.5, dpi = 150)

# The two bars corresponding to baseline b in the published Figure 9.
figure9_b$horizon_label <- factor(
  paste0(figure9_b$horizon, " quarters"),
  levels = c("4 quarters", "8 quarters")
)
p_bars <- ggplot(figure9_b, aes(horizon_label, estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  labs(
    title = "Figure 9: baseline specification b",
    subtitle = "Cumulated response of inflation to a temporary 1 ppt contractionary shock",
    x = NULL, y = "%"
  ) +
  theme_minimal(base_size = 11)
ggsave(outp("figure9_baseline_b_reviewed.png"), p_bars,
       width = 5.5, height = 4.5, dpi = 150)

# -----------------------------------------------------------------------------
# Four-variable IRF panel (order: cash rate, GDP, inflation, price level)
# -----------------------------------------------------------------------------
# The baseline VAR contains the raw policy SHOCK (not the cash rate), so the
# fourth panel is the monetary policy shock itself: a temporary innovation that
# reverts to zero. Its cumulation is the I(1) cash-rate-level proxy, which is
# permanent by construction - hence it is NOT comparable to the actual cash rate
# used in the LP and proxy-SVAR, and is deliberately excluded from the overlay.
var_labels <- c(
  mp_shock    = "Monetary policy shock (ppt)",
  real_gdp    = "Real GDP (% deviation)",
  inflation   = "Underlying inflation (quarterly ppt)",
  price_level = "Price level (% deviation)"
)
var_palette <- c(
  "Monetary policy shock (ppt)"          = "#4E6E58",
  "Real GDP (% deviation)"               = "#2C6E9B",
  "Underlying inflation (quarterly ppt)" = "#C1553B",
  "Price level (% deviation)"            = "#B23A48"
)

tidy_central <- function(x, response, varname) {
  data.frame(
    variable = varname,
    horizon  = 0:(nrow(x$irf$mpshock) - 1),
    estimate = x$irf$mpshock[, response],
    lower    = x$Lower$mpshock[, response],
    upper    = x$Upper$mpshock[, response]
  )
}

# Tidy central IRFs (also written to CSV for the combined overlay script).
var_irfs <- rbind(
  tidy_central(irf_shock,     "mpshock", "mp_shock"),
  tidy_central(irf_gdp_level, "ln_rgdp", "real_gdp"),
  tidy_central(irf_infl,      "infl",    "inflation"),
  tidy_central(irf_price,     "infl",    "price_level")
)
write.csv(var_irfs, outp("baseline_var_irfs.csv"), row.names = FALSE)

Hmax <- 20
var_plot_df <- var_irfs[var_irfs$horizon <= Hmax, ]
var_plot_df$variable_label <- factor(var_labels[var_plot_df$variable],
                                     levels = unname(var_labels))
var_markers <- subset(var_plot_df, horizon %in% c(4, 8))

theme_irf <- theme_minimal(base_size = 11) +
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

p_var_all <- ggplot(var_plot_df,
                    aes(horizon, estimate, colour = variable_label, fill = variable_label)) +
  geom_hline(yintercept = 0, colour = "#9aa0a6", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.0, lineend = "round") +
  geom_point(data = var_markers, shape = 21, colour = "white", size = 2.3, stroke = 0.7) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = var_palette, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = seq(0, Hmax, by = 4), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Baseline VAR: responses to a temporary +1 ppt policy shock",
    subtitle = "Recursive VAR(3), Cholesky order GDP-inflation-shock. Policy variable is the shock itself (temporary); price level = cumulated inflation.",
    x = "Quarters after shock", y = NULL,
    caption = "Bands are vars-package 90% percentile bootstrap intervals. Dots mark h = 4 and 8."
  ) +
  theme_irf
ggsave(outp("extension_baseline_var_all.png"), p_var_all, width = 10, height = 7, dpi = 160)

# -----------------------------------------------------------------------------
# Paper-text specification (not the specification used by the public Figure 9 loop)
# -----------------------------------------------------------------------------
paper_df <- var_df
paper_df$cum_shock <- cumsum(paper_df$mpshock)
# The do-file also recentres the cumulative series to the average cash rate. With
# an intercept this additive constant does not affect the VAR slopes or IRFs.
paper_ts <- ts(
  paper_df[, c("ln_rgdp", "infl", "cum_shock")],
  start = c(1991, 2), frequency = 4
)
fit_paper_text <- VAR(paper_ts, p = 3, type = "const")
paper_scale <- irf(
  fit_paper_text, impulse = "cum_shock", response = "cum_shock",
  n.ahead = 1, ortho = TRUE, boot = FALSE
)$irf$cum_shock[1, "cum_shock"]
paper_irf <- irf(
  fit_paper_text, impulse = "cum_shock", response = "infl",
  n.ahead = h, ortho = TRUE, cumulative = TRUE,
  boot = TRUE, runs = 1000, ci = 0.90, seed = 123
)
paper_price <- data.frame(
  horizon = 0:h,
  estimate = paper_irf$irf$cum_shock[, "infl"] / paper_scale * 100,
  lower = paper_irf$Lower$cum_shock[, "infl"] / paper_scale * 100,
  upper = paper_irf$Upper$cum_shock[, "infl"] / paper_scale * 100
)
cat("\n=== Paper-prose cumulative-shock specification (comparison only) ===\n")
print(paper_price[paper_price$horizon %in% c(4, 8), ], digits = 6, row.names = FALSE)

saveRDS(
  list(
    fit_var = fit_var, impact_scale = impact_scale,
    irf_price = irf_price, irf_infl = irf_infl,
    irf_gdp_level = irf_gdp_level, irf_gdp_do_coirf = irf_gdp_do_coirf,
    irf_shock = irf_shock,
    price_path = price_path, figure9_b = figure9_b, var_irfs = var_irfs,
    fit_paper_text = fit_paper_text, paper_price = paper_price,
    var_df = var_df
  ),
  outp("task2_results_reviewed.rds")
)
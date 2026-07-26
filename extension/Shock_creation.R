# Bishop & Tulip (2017), Equation (9): exact public-replication specification
#
# Delta r_m = alpha + phi*piF_m + gamma*yF_m + lambda*piFrev_m
#             + theta*yFrev_m + rho*uNowcast_m + beta*r_{m-1} + epsilon_m

library(readxl)
library(sandwich)
library(lmtest)
library(car)
library(ggplot2)
library(patchwork)

#Set your working directory here. 

setwd("C:/Users/Ben/Documents/UNSW course/Presentation/Australian replication/extension")

# All generated files are written to an "output" subfolder.
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outp <- function(f) file.path(out_dir, f)

q <- as.data.frame(read_excel("data.xlsx", sheet = "STATA_q"))
q <- q[order(q$year, q$quarter), ]
q$qidx <- q$year + (q$quarter - 1) / 4

reg_vars <- c(
  "cr_change", "yeinff_2", "yegdpf_2", "yeinff_2r",
  "yegdpf_2r", "unef_0", "cr_prev"
)

# Stata replication:
# newey ... if qtr>=tq(1991q1), lag(1) force
eligible <- q$qidx >= 1991.00
d <- q[eligible & complete.cases(q[, reg_vars]), ]
stopifnot(nrow(d) == 100)
stopifnot(d$year[1] == 1991, d$quarter[1] == 2)
stopifnot(d$year[nrow(d)] == 2016, d$quarter[nrow(d)] == 1)

fit <- lm(
  cr_change ~ yeinff_2 + yegdpf_2 + yeinff_2r +
    yegdpf_2r + unef_0 + cr_prev,
  data = d
)

# Exact analogue of Stata newey ..., lag(1): Bartlett HAC, one lag,
# no prewhitening, with finite-sample adjustment.
nw_vcov <- NeweyWest(fit, lag = 1, prewhite = FALSE, adjust = TRUE)
nw_test <- coeftest(fit, vcov. = nw_vcov)

cat("\n=== Equation (9): OLS coefficients and Newey-West HAC SEs ===\n")
print(nw_test)

# Stata's post-newey test is reported as a Wald chi-squared test. Explicitly
# request Chisq here; car::linearHypothesis otherwise defaults to an F test for lm.
joint_inf <- linearHypothesis(
  fit, c("yeinff_2 = 0", "yeinff_2r = 0"),
  vcov. = nw_vcov, test = "Chisq"
)
joint_gdp <- linearHypothesis(
  fit, c("yegdpf_2 = 0", "yegdpf_2r = 0"),
  vcov. = nw_vcov, test = "Chisq"
)
inf_sum_test <- linearHypothesis(
  fit, "yeinff_2 + yeinff_2r = 0",
  vcov. = nw_vcov, test = "Chisq"
)
gdp_sum_test <- linearHypothesis(
  fit, "yegdpf_2 + yegdpf_2r = 0",
  vcov. = nw_vcov, test = "Chisq"
)

cat("\n=== Joint Wald tests ===\n")
cat("Inflation forecast level and revision jointly zero:\n")
print(joint_inf)
cat("GDP forecast level and revision jointly zero:\n")
print(joint_gdp)

inf_sum <- coef(fit)["yeinff_2"] + coef(fit)["yeinff_2r"]
gdp_sum <- coef(fit)["yegdpf_2"] + coef(fit)["yegdpf_2r"]
cat(sprintf("\nInflation-coefficient sum = %.6f (paper: 0.19)\n", inf_sum))
cat(sprintf("GDP-coefficient sum       = %.6f (paper: 0.06)\n", gdp_sum))

r2 <- summary(fit)$r.squared
adj_r2 <- summary(fit)$adj.r.squared
dw <- dwtest(fit)
cat(sprintf("\nUnadjusted R-squared = %.6f\n", r2))
cat(sprintf("Adjusted R-squared   = %.6f\n", adj_r2))
cat(sprintf("Durbin-Watson        = %.6f\n", unname(dw$statistic)))
cat("Note: Table B1 labels 0.20 as adjusted R-squared, but 0.20 matches\n")
cat("the unadjusted R-squared in the released data; the adjusted value is about 0.153.\n")

# Table B1 audit
paper_targets <- data.frame(
  term = c("(Intercept)", "yeinff_2", "yegdpf_2", "yeinff_2r",
           "yegdpf_2r", "unef_0", "cr_prev"),
  paper_coef = c(-0.10, -0.05, 0.05, 0.24, 0.01, 0.00, 0.01),
  paper_se = c(0.19, 0.05, 0.03, 0.08, 0.05, 0.01, 0.02),
  stringsAsFactors = FALSE
)
results <- data.frame(
  term = rownames(nw_test),
  estimate = nw_test[, 1],
  hac_se = nw_test[, 2],
  stringsAsFactors = FALSE
)
audit <- merge(paper_targets, results, by = "term", all.x = TRUE, sort = FALSE)
audit <- audit[match(paper_targets$term, audit$term), ]
audit$coef_difference <- audit$estimate - audit$paper_coef
audit$se_difference <- audit$hac_se - audit$paper_se
cat("\n=== Audit against Table B1 ===\n")
print(audit, digits = 4, row.names = FALSE)

# The policy-shock series is the OLS residual. HAC affects inference only; it does
# not change fitted values or residuals.
shocks <- data.frame(
  year = d$year,
  quarter = d$quarter,
  date_inf_fore = d$date_inf_fore,
  shock = residuals(fit)
)
write.csv(shocks, outp("task1_shocks.csv"), row.names = FALSE)
cat(sprintf("\nSaved %d shocks: %dQ%d to %dQ%d\n",
            nrow(shocks), shocks$year[1], shocks$quarter[1],
            shocks$year[nrow(shocks)], shocks$quarter[nrow(shocks)]))

# Optional Figure 4-style diagnostic. Handle either Date/POSIX or Excel serial dates.
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(x, origin = "1899-12-30")
}
plot_df <- data.frame(
  date = to_date(d$date_inf_fore),
  actual = d$cr_change,
  fitted = fitted(fit),
  residual = residuals(fit)
)

p_top <- ggplot(plot_df, aes(date)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(y = actual, colour = "Actual"), linewidth = 0.7) +
  geom_line(aes(y = fitted, colour = "Fitted"), linewidth = 0.7) +
  labs(title = "Changes in the cash rate", x = NULL, y = "ppt", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

p_bottom <- ggplot(plot_df, aes(date, residual)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 55) +
  labs(title = "Equation (9) residuals", x = NULL, y = "ppt") +
  theme_minimal(base_size = 11)

ggsave(outp("extension_monetary_shocks.png"), p_top / p_bottom,
       width = 8.5, height = 6.5, dpi = 150)

saveRDS(
  list(fit = fit, nw_vcov = nw_vcov, nw_test = nw_test,
       joint_inf = joint_inf, joint_gdp = joint_gdp,
       inf_sum_test = inf_sum_test, gdp_sum_test = gdp_sum_test,
       r2 = r2, adj_r2 = adj_r2, dw = dw, audit = audit,
       shocks = shocks),
  outp("task1_results_reviewed.rds")
)
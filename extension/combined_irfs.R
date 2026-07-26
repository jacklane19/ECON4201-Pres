# =============================================================================
# Combined overlay of CENTRAL IRFs across the three identifications:
#   1. Baseline VAR        -> baseline_var_irfs.csv        (from Baseline_var.R)
#   2. Direct-shock LP      -> direct_shock_lp_results.csv  (from direct_shock_lp.R)
#   3. Proxy SVAR           -> proxy_svar_irfs.csv          (from the proxy-SVAR script)
#
# Run the three specification scripts first so those CSVs exist, then run this.
# Panels are ordered cash rate, GDP, inflation, price level; one line per spec.
# =============================================================================

library(dplyr)
library(ggplot2)

H_overlay <- 16   # common horizon across the three specifications

# Inputs and outputs live in the "output" subfolder.
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
outp <- function(f) file.path(out_dir, f)

# Each spec stores its variable key under a different column name. Look in the
# output folder first, then fall back to the working directory.
read_spec <- function(path, var_col, spec_name) {
  resolved <- if (file.exists(outp(path))) outp(path) else path
  if (!file.exists(resolved)) {
    stop("Missing '", path, "' (looked in '", out_dir,
         "/' and the working directory). Run the corresponding specification script first.")
  }
  df <- read.csv(resolved, stringsAsFactors = FALSE)
  data.frame(
    spec     = spec_name,
    variable = df[[var_col]],
    horizon  = df$horizon,
    estimate = df$estimate,
    stringsAsFactors = FALSE
  )
}

combined <- bind_rows(
  read_spec("baseline_var_irfs.csv",       "variable", "Baseline VAR"),
  read_spec("direct_shock_lp_results.csv", "outcome",  "Direct-shock LP"),
  read_spec("proxy_svar_irfs.csv",         "response", "Proxy SVAR")
) %>%
  filter(
    # Cash rate is deliberately excluded. The baseline VAR has no cash-rate
    # variable, so its "cash rate" would be the cumulated-shock (I(1)) proxy,
    # which is permanent by construction and NOT comparable to the actual,
    # mean-reverting cash rate used in the LP and proxy-SVAR. Only responses that
    # are like-for-like across the three identifications are overlaid.
    variable %in% c("real_gdp", "inflation", "price_level"),
    horizon >= 0, horizon <= H_overlay
  )

# Panel order and labels (GDP, inflation, price level).
var_labels <- c(
  real_gdp    = "Real GDP (% deviation)",
  inflation   = "Underlying inflation (quarterly ppt)",
  price_level = "Price level (% deviation)"
)
combined$variable_label <- factor(var_labels[combined$variable],
                                  levels = unname(var_labels))
combined$spec <- factor(combined$spec,
                        levels = c("Baseline VAR", "Direct-shock LP", "Proxy SVAR"))

spec_palette <- c(
  "Baseline VAR"    = "#2C6E9B",
  "Direct-shock LP" = "#C1553B",
  "Proxy SVAR"      = "#4E6E58"
)

theme_overlay <- theme_minimal(base_size = 11) +
  theme(
    plot.title           = element_text(face = "bold", size = 14, colour = "#16181c"),
    plot.subtitle        = element_text(size = 9.5, colour = "#5f6368", margin = margin(b = 6)),
    plot.caption         = element_text(size = 8, colour = "#9aa0a6", hjust = 0),
    strip.text           = element_text(face = "bold", size = 11, colour = "#20232a", hjust = 0),
    axis.title           = element_text(size = 8.5, colour = "#5f6368"),
    axis.text            = element_text(size = 8.5, colour = "#5f6368"),
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_line(colour = "#eef0f2", linewidth = 0.5),
    axis.ticks           = element_blank(),
    legend.position      = "top",
    legend.title         = element_blank(),
    legend.justification = "left",
    panel.spacing        = unit(1.1, "lines"),
    plot.margin          = margin(8, 14, 8, 8)
  )

p_combined <- ggplot(combined, aes(horizon, estimate, colour = spec)) +
  geom_hline(yintercept = 0, colour = "#9aa0a6", linewidth = 0.35) +
  geom_line(linewidth = 1.05, lineend = "round") +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = spec_palette) +
  scale_x_continuous(breaks = seq(0, H_overlay, by = 4),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "Central IRFs across three identifications",
    subtitle = "Response to a +1 ppt monetary policy shock; point estimates only (no confidence bands).",
    x = "Quarters after shock", y = NULL,
    caption = paste0(
      "Baseline VAR, direct-shock LP and proxy-SVAR, horizons 0-", H_overlay,
      ". Cash rate omitted: not comparable across identifications (the VAR has no cash-rate variable)."
    )
  ) +
  theme_overlay

ggsave(outp("combined_irfs.png"), p_combined, width = 12, height = 4.6, dpi = 170, bg = "white")
write.csv(combined, outp("combined_irfs_central.csv"), row.names = FALSE)

cat("Saved", outp("combined_irfs.png"), "and", outp("combined_irfs_central.csv"), "\n")

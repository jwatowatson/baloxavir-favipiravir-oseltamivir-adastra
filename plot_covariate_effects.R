library(dplyr)
library(ggplot2)

main_mod <- 1
load(ff[main_mod])  # loads `out`, `stan_inputs`, etc.

# Helper: summarize draws to 2.5%, 50%, 97.5%
qsum <- function(df, probs = c(0.025, 0.5, 0.975)) {
  as.data.frame(t(apply(df, 2, stats::quantile, probs = probs))) |>
    `colnames<-`(c("Low", "Med", "Up")) |>
    tibble::rownames_to_column("cov")
}

# 1) Extract posterior draws ---------------------------------------------------------
# Slope
slope_draws <- rstan::extract(out, "slope_coefs")$slope_coefs |> as.data.frame()
colnames(slope_draws) <- colnames(stan_inputs[[1]]$cov_matrices$X_slope[[1]])
slope_summary <- qsum(slope_draws) |> mutate(type = "Slope")
slope_summary
# Intercept
int_draws <- rstan::extract(out, "intercept_coefs")$intercept_coefs |> as.data.frame()
colnames(int_draws) <- colnames(stan_inputs[[1]]$cov_matrices$X_int[[1]])
int_summary <- qsum(int_draws) |> mutate(type = "Intercept")
int_summary
# 2) Back-transform age effect (SD-wise -> per +10 years) ----------------------------
sd_age <- adastra_dat_analysis$SD_age[1]
int_summary <- int_summary |>
  mutate(across(c(Low, Med, Up), ~ ifelse(cov == "Age_scaled", . * 10 / sd_age, .)))

# 3) Transform to interpretable scales ----------------------------------------------
# - Slope: exp(coef) then format as percent change (via `formatter`)
# - Intercept: 10^(coef) then format as percent change (via `formatter`)
slope_summary <- slope_summary |> mutate(across(c(Low, Med, Up), ~ formatter(exp(.))))
int_summary  # <- int_summary   |> mutate(across(c(Low, Med, Up), ~ formatter(10^(.))))

# Bind summaries
coefs_summary <- bind_rows(slope_summary, int_summary) |>
  mutate(cov = factor(cov))

# 4) Optional: relabel covariates for nicer display ---------------------------------
label_map <- c(
  Age_scaled = "Age (+10 years)",
  vaccineYes  = "Vaccine status (yes vs no)",
  symptomDay = "Symptom day (+1)",
  Study_time = "Study time (study-centered)",
  fluTypeB   = "Influenza B (vs A)",
  Sex = "Sex (male vs female)",
  Sitela008 = "Site: Laos (vs Brazil)",
  Sitenp003 = "Site: Nepal (vs Brazil)",
  Siteth001 = "Site: Thailand (vs Brazil)"
)
coefs_summary <- coefs_summary |>
  mutate(cov = forcats::fct_relabel(cov, ~ label_map[.x] %||% .x))

# 5) Plotting ------------------------------------------------------------------------
base_theme <- theme_bw(base_size = 14) +
  theme(
    strip.text   = element_text( face = "bold"),
    axis.title   = element_text(face = "bold"),
    plot.title   = element_blank(),
    axis.text    = element_text(),
    legend.position = "bottom",
    panel.spacing   = grid::unit(1, "lines")
  )

G_cov1 <- coefs_summary |> filter(type == "Slope") |>
  ggplot(aes(x = cov, y = Med)) +
  geom_errorbar(aes(ymin = Low, ymax = Up), width = 0, alpha = 0.75, linewidth = 1.5, col = "#1D2B53") +
  geom_point(size = 3.5, col = "#1D2B53") +
  geom_hline(yintercept = 0, linetype = "32", col = "#CF0F0F") +
  coord_flip() +
  labs(title = "B) Effects on viral clearance rate",
       x = "", y = "Change in viral clearance rate (%)") +
  base_theme

G_cov2 <- coefs_summary |> filter(type == "Intercept") |>
  ggplot(aes(x = cov, y = Med)) +
  geom_errorbar(aes(ymin = Low, ymax = Up), width = 0, alpha = 0.75, linewidth = 1.5, col = "#1D2B53") +
  geom_point(size = 3.5, col = "#1D2B53") +
  geom_hline(yintercept = 0, linetype = "32", col = "#CF0F0F") +
  coord_flip() +
  labs(title = "A) Effects on baseline viral density",
       x = "", y = "Change in admission viral densities (log10 copies/mL)") +
  base_theme

G_cov1
G_cov2

G_combined <- ggarrange(G_cov2, G_cov1,
          ncol = 1, nrow = 2,
          labels = c("", ""),
          align = "v")
G_combined

plot_name <- here("Plots", intervention, paste0("mITT_", mITT_threshold, "_covariate_effects.png"))
png(plot_name, width = 8, height = 8, units = "in", res = 350)
print(G_combined)
dev.off()

# Load packages ####

library(tidyverse)
library(amt)

# Load models ####

issa_spring <- readRDS("output/iSSA_spring_2024-06-19.rds")
issa_fall <- readRDS("output/iSSA_fall_2024-06-19.rds")

# Inverse variance-weighted regression ####

# Format model output
issa_spring <- issa_spring %>%
  mutate(coef = map(issa, function(x) broom::tidy(x$model)))

issa_fall <- issa_fall %>%
  mutate(coef = map(issa, function(x) broom::tidy(x$model)))

# Unnest coefficients
coefs_spring <- issa_spring %>%
  select(-issa) %>%
  unnest(cols = coef)

coefs_fall <- issa_fall %>%
  select(-issa) %>%
  unnest(cols = coef)

# Nest data frame so that I have one per term (covariate)
coefs_iv_spring <- coefs_spring %>%
  nest(data = -term)

coefs_iv_fall <- coefs_fall %>%
  nest(data = -term)

# Get rid of NAs (parameters that could not be estimated, e.g., individuals
# whose steps never fell in a particular land cover)
coefs_iv_spring_clean <- coefs_iv_spring %>%
  mutate(data = lapply(data, function(x) {
    x %>%
      filter(!is.na(estimate) & !is.nan(statistic))
  }))

coefs_iv_fall_clean <- coefs_iv_fall %>%
  mutate(data = lapply(data, function(x) {
    x %>%
      filter(!is.na(estimate) & !is.nan(statistic))
  }))

# Fit linear model and get predictions
coefs_iv_spring_mod <- coefs_iv_spring_clean %>%
  mutate(iv = lapply(data, function(x) {
    mod <- lm(estimate ~ 1, data = x, weights = 1/(std.error)^2)
    return(mod)
  })) %>%
  mutate(pred = lapply(iv, function(x) {
    pred <- predict(x,
                    newdata = data.frame(dummy = 1),
                    se.fit = TRUE)
    est <- data.frame(mean = pred$fit,
                      lwr = pred$fit - 1.96 * pred$se.fit,
                      upr = pred$fit + 1.96 * pred$se.fit)
    return(est)
  }))

coefs_iv_fall_mod <- coefs_iv_fall_clean %>%
  mutate(iv = lapply(data, function(x) {
    mod <- lm(estimate ~ 1, data = x, weights = 1/(std.error)^2)
    return(mod)
  })) %>%
  mutate(pred = lapply(iv, function(x) {
    pred <- predict(x,
                    newdata = data.frame(dummy = 1),
                    se.fit = TRUE)
    est <- data.frame(mean = pred$fit,
                      lwr = pred$fit - 1.96 * pred$se.fit,
                      upr = pred$fit + 1.96 * pred$se.fit)
    return(est)
  }))

# Unnest

iv_spring <- coefs_iv_spring_mod %>%
  dplyr::select(term, pred) %>%
  unnest(cols = pred)

iv_fall <- coefs_iv_fall_mod %>%
  dplyr::select(term, pred) %>%
  unnest(cols = pred)

# Plot parameter estimates ####

term_pretty_spring <- c("Distance to summer range : Cos(turning angle)",
                 "Distance to summer range",
                 "NDVI gain : Starting NDVI",
                 "NDVI gain",
                 "Distance to roads",
                 "Elevation gain (quadratic) : Starting elevation",
                 "Elevation gain (linear) : Starting elevation",
                 "Elevation gain (quadratic)",
                 "Elevation gain (linear)",
                 "Bare",
                 "Shrubland",
                 "Grassland",
                 "Developed",
                 "Agricultural")

iv_spring %>%
  mutate(term_order = c("T. Agricultural",
                        "S. Developed",
                        "R. Grassland",
                        "Q. Shrubland",
                        "P. Bare",
                        "O. Elevation gain (linear)",
                        "N. Elevation gain (quadratic)",
                        "I. Distance to roads",
                        "F. Distance to summer range",
                        "H. NDVI gain",
                        "D. Step length",
                        "C. Log(step length)",
                        "B. Cos(turning angle)",
                        "M. Elevation gain (linear) : Starting elevation",
                        "L. Elevation gain (quadratic) : Starting elevation",
                        "E. Distance to summer range : Cos(turning angle)",
                        "G. NDVI gain : Starting NDVI",
                        "A. Log(step length : Cos(turning angle)")) %>%
  filter(!term %in% c("cos_ta_", "log_sl_", "sl_", "log_sl_:cos_ta_")) %>%
  ggplot(aes(x = term_order, y = mean, color = term_order)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = " ", y = "log-RSS") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(size = 12)) +
  coord_flip() +
  scale_x_discrete(labels = term_pretty_spring) +
  # scale_color_manual(name = "Covariate",
  #                    values = c(rep("#F8AF67", 2), ##EC8447
  #                               rep("#767171", 9))) +
  scale_color_manual(name = "Covariate",
                     values = rep("#767171", 14))

ggsave("output/IVWR_log-RSS_spring.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)
ggsave("output/IVWR_log-RSS_spring_land-cover-highlight.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)
ggsave("output/IVWR_log-RSS_spring_elevation-highlight.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)
ggsave("output/IVWR_log-RSS_spring_roads-highlight.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)
ggsave("output/IVWR_log-RSS_spring_NDVI-highlight.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)
ggsave("output/IVWR_log-RSS_spring_dist-to-range-highlight.tiff",
       dpi = 300, compression = "lzw", width = 8, height = 3)

term_pretty_fall <- c("Distance to winter range : Cos(turning angle)",
                        "Distance to winter range",
                        "NDVI gain : Starting NDVI",
                        "NDVI gain",
                        "Distance to roads",
                        "Elevation gain (quadratic) : Starting elevation",
                        "Elevation gain (linear) : Starting elevation",
                        "Elevation gain (quadratic)",
                        "Elevation gain (linear)",
                        "Bare",
                        "Shrubland",
                        "Grassland",
                        "Developed",
                        "Agricultural")

iv_fall %>%
  mutate(term_order = c("T. Agricultural",
                        "S. Developed",
                        "R. Grassland",
                        "Q. Shrubland",
                        "P. Bare",
                        "O. Elevation gain (linear)",
                        "N. Elevation gain (quadratic)",
                        "I. Distance to roads",
                        "F. Distance to winter range",
                        "H. NDVI gain",
                        "D. Step length",
                        "C. Log(step length)",
                        "B. Cos(turning angle)",
                        "M. Elevation gain (linear) : Starting elevation",
                        "L. Elevation gain (quadratic) : Starting elevation",
                        "E. Distance to winter range : Cos(turning angle)",
                        "G. NDVI gain : Starting NDVI",
                        "A. Log(step length : Cos(turning angle)")) %>%
  filter(!term %in% c("cos_ta_", "log_sl_", "sl_", "log_sl_:cos_ta_")) %>%
  ggplot(aes(x = term_order, y = mean, color = term_order)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = " ", y = "log-RSS") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(size = 12)) +
  coord_flip() +
  scale_x_discrete(labels = term_pretty_fall) +
  # scale_color_manual(name = "Covariate",
  #                    values = c(rep("#F8AF67", 2), ##EC8447
  #                               rep("#767171", 9))) +
  scale_color_manual(name = "Covariate",
                     values = rep("#767171", 14))

library(tidyverse)
library(tseries)
library(urca)
library(vars)
library(lmtest)
library(sandwich)

# ============================================================================
# 1. SETUP
# ============================================================================

# Clean workspace (optional for scripts; omit in Rmd)
rm(list = ls())

# ============================================================================
# 2. DATA LOAD (LONG FORMAT)
# ============================================================================

df <- read_csv('data/series.csv') %>% 
  arrange(unique_id, ds)

# ============================================================================
# 3. EXPLORATORY PLOTS
# ============================================================================

plot_facets <- function(df_long, value_col, ylabel, title) {
  ggplot(df_long, aes(x = ds, y = .data[[value_col]])) +
    geom_line() +
    facet_wrap(~unique_id) +
    labs(title = title, x = NULL, y = ylabel) +
    theme_minimal(base_size = 11)
}

plot_facets(df, "y", "Indice", "Indice base 2018")

df <- df %>%
  group_by(unique_id) %>%
  mutate(yoy = (y / lag(y, 12) - 1) * 100) %>%
  ungroup()

plot_facets(df, "yoy", "Cambio porcentual anual (%)", "Variacion anual (%)")

# ============================================================================
# 4. STATIONARITY TESTS
# ============================================================================

adf_test <- function(series) {
  test <- adf.test(na.omit(series))
  tibble(
    adf_stat = unname(test$statistic),
    p_value = test$p.value,
    n_lags = unname(test$parameter)
  )
}

kpss_test <- function(series, null = "Trend") {
  test <- kpss.test(na.omit(series), null = null)
  tibble(
    kpss_stat = unname(test$statistic),
    kpss_p = test$p.value,
    kpss_lags = unname(test$parameter)
  )
}

adf_df <- df %>%
  group_by(unique_id) %>%
  group_modify(~adf_test(.x$y)) %>%
  ungroup()

kpss_df <- df %>%
  group_by(unique_id) %>%
  group_modify(~kpss_test(.x$y, null = "Trend")) %>%
  ungroup()

adf_df
kpss_df

# ============================================================================
# 5. DIFFERENCES
# ============================================================================

df <- df %>%
  group_by(unique_id) %>%
  mutate(dy = y - lag(y)) %>%
  ungroup()

plot_facets(df, "dy", "Delta y", "Delta y, indices diferenciados")

adf_dy <- df %>%
  group_by(unique_id) %>%
  group_modify(~adf_test(.x$dy)) %>%
  ungroup()

kpss_dy <- df %>%
  group_by(unique_id) %>%
  group_modify(~kpss_test(.x$dy, null = "Level")) %>%
  ungroup()

adf_dy
kpss_dy

# ============================================================================
# 6. VAR AND COINTEGRATION
# ============================================================================

df_wide <- df %>%
  dplyr::select(ds, unique_id, y) %>%
  tidyr::pivot_wider(names_from = unique_id, values_from = y) %>%
  dplyr::arrange(ds) %>%
  tidyr::drop_na(igae, consumo, inversion)

start_year <- lubridate::year(min(df_wide$ds))
start_month <- lubridate::month(min(df_wide$ds))

Y_ts <- ts(
  df_wide %>% dplyr::select(igae, consumo, inversion),
  start = c(start_year, start_month),
  frequency = 12
)

lag_sel <- VARselect(Y_ts, lag.max = 12, type = "const")
lag_sel$selection

K <- as.integer(lag_sel$selection["AIC(n)"])

joh_trace <- ca.jo(Y_ts, type = "trace", ecdet = "const", K = K)
summary(joh_trace)

joh_eigen <- ca.jo(Y_ts, type = "eigen", ecdet = "const", K = K)
summary(joh_eigen)

# ============================================================================
# 7. LOG TRANSFORM AND REGRESSION
# ============================================================================

df <- df %>%
  mutate(ly = log(y))

plot_facets(df, "ly", "Log", "Series en log")

df_wide_ly <- df %>%
  dplyr::filter(is.finite(ly)) %>%
  dplyr::select(ds, unique_id, ly) %>%
  tidyr::pivot_wider(names_from = unique_id, values_from = ly) %>%
  dplyr::arrange(ds) %>%
  tidyr::drop_na(igae, consumo, inversion)

model <- lm(igae ~ consumo + inversion, data = df_wide_ly)
coeftest(model, vcov = NeweyWest(model, lag = 12, prewhite = FALSE))

summary(model)

# ==============================================================================
# PROJECT: Finnish ENergy Market Day-Ahead Electricity Price Forecasting
# DESCRIPTION: Two-stage SARIMAX + eGARCH-t model with Fourier terms for 
#              double seasonality (daily/weekly) and volatility clustering.
# AUTHOR: Enio Hoxha
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
# Install required packages if not present
required_packages <- c("tidyverse", "lubridate", "forecast", "rugarch", 
                       "tseries", "lmtest", "gridExtra", "FinTS")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(forecast)
  library(rugarch)
  library(tseries)
  library(lmtest)
  library(gridExtra)
  library(FinTS)
})
DATA_FILE <- "Finnish_Energy_Market_Weather_Data.csv"
KFAS_FILE <- "kfas_benchmark_forecast.csv" # Optional: KFAS State-Space benchmark

if(!file.exists(DATA_FILE)) {
  stop("CRITICAL ERROR: Dataset not found. Please ensure 'Finnish_Energy_Market_Weather_Data.csv' is in the working directory.")
}

# ------------------------------------------------------------------------------
# 2. DATA LOADING & FEATURE ENGINEERING
# ------------------------------------------------------------------------------
cat("\n[1/6] Loading and Preprocessing Data\n")

df <- read_csv(DATA_FILE, show_col_types = FALSE) %>%
  mutate(datetime = as_datetime(datetime)) %>%
  arrange(datetime) %>%
  tidyr::fill(everything(), .direction = "down") %>%
  mutate(
    price = sp,
    # Standardize Exogenous Regressor (Wind Speed)
    wind_std = as.numeric(scale(wind_speed))
  ) %>%
  slice(-1) # Remove first row in case of lagged NAs

# Create Time Series Object with 168h (Weekly) Seasonality
msts_price <- msts(df$price, seasonal.periods = 168)

# Generate Fourier Terms (K=2) to capture weekly seasonality and remove the "wave"
weekly_fourier <- fourier(msts_price, K = 2)

# Combine active regressors
xreg_full <- cbind(wind_std = df$wind_std, weekly_fourier)

# ------------------------------------------------------------------------------
# 3. FULL SAMPLE ESTIMATION: SARIMAX + eGARCH
# ------------------------------------------------------------------------------
cat("\n[2/6] Estimating Mean Equation: SARIMAX(1,1,2)(0,1,1)[24] + Fourier\n")

fit_sarimax <- Arima(df$price, 
                     order = c(1, 1, 2), 
                     seasonal = list(order = c(0, 1, 1), period = 24),
                     xreg = xreg_full, 
                     method = "CSS-ML")

z_mean <- residuals(fit_sarimax)

cat("\n[3/6] Estimating Variance Equation: eGARCH(1,1)-Student-t\n")

spec_egarch <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std"
)

fit_egarch <- ugarchfit(spec = spec_egarch, data = z_mean, solver = "hybrid")

# ------------------------------------------------------------------------------
# 4. RIGOROUS ECONOMETRIC DIAGNOSTICS
# ------------------------------------------------------------------------------
cat("\n[4/6] Running Diagnostics\n")

std_resid <- residuals(fit_egarch, standardize = TRUE)
std_resid_sq <- std_resid^2 

# Statistical Tests
lb_mean <- Box.test(std_resid, lag = 48, type = "Ljung-Box")
lb_vol  <- Box.test(std_resid_sq, lag = 48, type = "Ljung-Box")
arch_lm <- ArchTest(as.numeric(std_resid), lags = 24)

cat("--- DIAGNOSTIC RESULTS ---\n")
cat(sprintf("Ljung-Box (Mean Dynamics) p-value:       %.4f\n", lb_mean$p.value))
cat(sprintf("Ljung-Box (Volatility Dynamics) p-value: %.4f\n", lb_vol$p.value))
cat(sprintf("ARCH-LM Test p-value:                    %.4f\n", arch_lm$p.value))

# Diagnostic Plots
p_acf_mean <- ggAcf(std_resid, lag.max = 168) + 
  ggtitle("ACF: Std Residuals (Weekly Wave Removed)") + theme_minimal()

p_acf_vol <- ggAcf(std_resid_sq, lag.max = 168) + 
  ggtitle("ACF: Squared Std Residuals (ARCH Effects Removed)") + theme_minimal()

df_shape <- coef(fit_egarch)["shape"]
qq_data <- data.frame(Sample = sort(as.numeric(std_resid)), 
                      Theoretical = qt(ppoints(length(std_resid)), df = df_shape))

p_qq <- ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1) +
  ggtitle(paste0("Q-Q Plot (Student-t, df=", round(df_shape, 2), ")")) + theme_minimal()

# Save diagnostic plot to disk (useful for GitHub README)
diag_plot <- grid.arrange(p_acf_mean, p_acf_vol, p_qq, nrow = 3)
ggsave("diagnostic_plots.png", diag_plot, width = 8, height = 10)


# ------------------------------------------------------------------------------
# 5. EXPANDING WINDOW BACKTEST (PSEUDO REAL-TIME)
# ------------------------------------------------------------------------------
cat("\n[5/6] Running Expanding Window Backtest (Last 7 Days)\n")

N_TEST_DAYS <- 7
HORIZON <- 24
total_obs <- nrow(df)
test_size <- N_TEST_DAYS * HORIZON
start_index <- total_obs - test_size

results_list <- list()

for (i in 0:(N_TEST_DAYS - 1)) {
  
  train_end_idx <- start_index + (i * HORIZON)
  train_indices <- 1:train_end_idx
  test_indices  <- (train_end_idx + 1):(train_end_idx + HORIZON)
  
  if(max(test_indices) > total_obs) break
  
  train_data <- df[train_indices, ]
  test_data  <- df[test_indices, ]
  
  # Dynamic Fourier Terms for current training window
  msts_train <- msts(train_data$price, seasonal.periods = 168)
  fourier_train <- fourier(msts_train, K = 2)
  fourier_test  <- fourier(msts_train, K = 2, h = HORIZON)
  
  xreg_train_loop <- cbind(wind_std = train_data$wind_std, fourier_train)
  xreg_test_loop  <- cbind(wind_std = test_data$wind_std, fourier_test)
  
  # Estimate SARIMAX
  fit_sarimax_loop <- tryCatch({
    Arima(train_data$price, order = c(1, 1, 2), seasonal = list(order = c(0, 1, 1), period = 24),
          xreg = xreg_train_loop, method = "CSS-ML")
  }, error = function(e) auto.arima(train_data$price, xreg=xreg_train_loop))
  
  pred_mean <- as.numeric(forecast(fit_sarimax_loop, xreg = xreg_test_loop, h = HORIZON)$mean)
  
  # Estimate eGARCH
  fit_garch_loop <- tryCatch({
    ugarchfit(spec = spec_egarch, data = residuals(fit_sarimax_loop), solver = 'hybrid')
  }, error = function(e) NULL)
  
  if(!is.null(fit_garch_loop)) {
    pred_sigma <- as.numeric(sigma(ugarchforecast(fit_garch_loop, n.ahead = HORIZON)))
    t_crit <- qt(0.975, df = coef(fit_garch_loop)["shape"])
  } else {
    pred_sigma <- rep(sd(residuals(fit_sarimax_loop)), HORIZON)
    t_crit <- qnorm(0.975)
  }
  
  results_list[[i + 1]] <- data.frame(
    datetime = test_data$datetime,
    Actual = test_data$price,
    Forecast = pred_mean,
    Lower_95 = pred_mean - t_crit * pred_sigma,
    Upper_95 = pred_mean + t_crit * pred_sigma
  )
}

backtest_results <- bind_rows(results_list)
rmse_val <- sqrt(mean((backtest_results$Actual - backtest_results$Forecast)^2))
mae_val <- mean(abs(backtest_results$Actual - backtest_results$Forecast))
cat(sprintf("Backtest Results -> RMSE: %.2f | MAE: %.2f\n", rmse_val, mae_val))


# ------------------------------------------------------------------------------
# 6. TRUE OUT-OF-SAMPLE FORECAST & MODEL COMPARISON (TARGET: MARCH 17)
# ------------------------------------------------------------------------------
cat("\n[6/6] Generating Final Out-of-Sample Forecast and Comparison Plot\n")

# Prepare Future Regressors (Assuming Persistence for Wind)
last_date <- as.Date(tail(df$datetime, 1))
weather_source <- df %>% filter(date(datetime) == last_date) %>% head(24)

# Future Fourier Terms
future_fourier <- fourier(msts_price, K = 2, h = 24)
future_xreg <- cbind(wind_std = weather_source$wind_std, future_fourier)

# Point Forecast
fcast_final <- forecast(fit_sarimax, xreg = future_xreg, h = 24)
pred_mean_final <- as.numeric(fcast_final$mean)

# Volatility Forecast
fcast_egarch_final <- ugarchforecast(fit_egarch, n.ahead = 24)
pred_sigma_final <- as.numeric(sigma(fcast_egarch_final))
t_crit_final <- qt(0.975, df = coef(fit_egarch)["shape"])
t_crit_80_final <- qt(0.90, df = coef(fit_egarch)["shape"])

future_dates <- seq(from = tail(df$datetime, 1) + hours(1), by = "1 hour", length.out = 24)

forecast_df <- data.frame(
  datetime = future_dates,
  Forecast_Price = pred_mean_final,
  Lower_95 = pred_mean_final - t_crit_final * pred_sigma_final,
  Upper_95 = pred_mean_final + t_crit_final * pred_sigma_final,
  Lower_80 = pred_mean_final - t_crit_80_final * pred_sigma_final,
  Upper_80 = pred_mean_final + t_crit_80_final * pred_sigma_final,
  Type = "SARIMAX+eGARCH"
)

# Integrate KFAS State-Space Data
if(file.exists(KFAS_FILE)) {
  df_ssm <- read_csv(KFAS_FILE, show_col_types = FALSE)
  num_cols <- select(df_ssm, where(is.numeric))
  if(ncol(num_cols) > 0 && nrow(num_cols) == 24) {
    ssm_df <- data.frame(
      datetime = future_dates, Forecast_Price = num_cols[[1]],
      Lower_95 = NA, Upper_95 = NA, Lower_80 = NA, Upper_80 = NA,
      Type = "KFAS State-Space"
    )
    forecast_df <- bind_rows(forecast_df, ssm_df)
  }
}

# Prepare Plot Data
history_plot <- df %>% select(datetime, price) %>%
  mutate(Forecast_Price = price, Lower_95 = NA, Upper_95 = NA, 
         Lower_80 = NA, Upper_80 = NA, Type = "Observed Price") %>%
  rename(Actual_Price = price) %>% tail(24 * 4)

plot_data <- bind_rows(history_plot, forecast_df)

plot_colors <- c("Observed Price" = "black", "SARIMAX+eGARCH" = "blue", "KFAS State-Space" = "forestgreen")
plot_linetypes <- c("Observed Price" = "solid", "SARIMAX+eGARCH" = "dashed", "KFAS State-Space" = "solid")

# Generate Professional Plot
p_final <- ggplot(plot_data, aes(x = datetime)) +
  geom_ribbon(data = subset(plot_data, Type == "SARIMAX+eGARCH"),
              aes(ymin = Lower_95, ymax = Upper_95, fill = "95% Prediction Interval"), alpha = 0.15) +
  geom_ribbon(data = subset(plot_data, Type == "SARIMAX+eGARCH"),
              aes(ymin = Lower_80, ymax = Upper_80, fill = "80% Prediction Interval"), alpha = 0.25) +
  geom_line(aes(y = Forecast_Price, color = Type, linetype = Type), size = 1) +
  geom_vline(xintercept = as.numeric(tail(df$datetime, 1)), linetype = "dotted", color = "gray30") +
  scale_color_manual(name = "Model / Series", values = plot_colors) +
  scale_linetype_manual(name = "Model / Series", values = plot_linetypes) +
  scale_fill_manual(name = "Uncertainty (SARIMAX)", values = c("95% Prediction Interval" = "firebrick", "80% Prediction Interval" = "red")) +
  labs(
    title = "Day-Ahead Electricity Price Forecast Comparison",
    subtitle = "Econometric (SARIMAX+eGARCH) vs State-Space (KFAS) Approaches",
    y = "Price (EUR/MWh)", x = "Date/Time (UTC)",
    caption = "Confidence intervals generated via Student-t distribution."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"), legend.box = "vertical")

print(p_final)
ggsave("forecast_comparison.png", p_final, width = 10, height = 6)

cat("\nPipeline Execution Complete! Visualizations saved to disk.\n")
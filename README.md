# finnish-electricity-price-forecasting
Econometric time-series forecasting of hourly day-ahead electricity prices in the Finnish market (Nord Pool) using SARIMAX-eGARCH and State Space models in R.
# Day-Ahead Electricity Price Forecasting: SARIMAX-eGARCH vs. State Space Models

[![Language](https://img.shields.io/badge/Language-R-276DC3.svg?logo=r)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-Completed-success)]()
[![Domain](https://img.shields.io/badge/Domain-Energy%20Trading%20%7C%20Nord%20Pool-blue)]()
[![Methodology](https://img.shields.io/badge/Methodology-Time%20Series%20Econometrics-orange)]()

An econometric time-series framework for forecasting hourly day-ahead electricity prices in the **Finnish Energy Market (Nord Pool)**, accounting for double seasonality, price spikes, and conditional volatility.

The repository includes a fully reproducible econometric pipeline in R comparing a **two-stage hybrid SARIMAX-eGARCH model** against a dynamic **State Space Model (Kalman Filter)** benchmark.

**Project Slides:** [Download Presentation PDF](docs/finnish-electricity-price-forecasting.pdf)

---

## Project Architecture & Stylized Facts

Wholesale electricity prices feature complex dynamics:
* **High-frequency double seasonality:** Interlocking daily (24h) and weekly (168h) cycles.
* **Extreme price spikes & leptokurtosis:** Fat-tailed distributions driven by supply-demand inelasticity.
* **Volatility clustering (ARCH effects):** Periods of market turbulence alternating with relative calm.

### Modeled Architectures
1. **Mean Equation (SARIMAX):** `SARIMA(1,1,2)(0,1,1)[24]` augmented with harmonic **Fourier terms ($K=2$)** to capture long-cycle weekly patterns, alongside exogenous standardized wind generation regressors ($p < 0.001$).
2. **Variance Equation (eGARCH):** `eGARCH(1,1)` under a **Student-t distribution** to capture asymmetric leverage effects and fat-tailed shocks.
3. **Benchmark (State Space Model / KFAS):** Time-varying local trend and stochastic seasonal components estimated via the Kalman Filter.

---

## Econometric Methodology

### 1. Dual-Seasonality Resolution
Standard Box-Jenkins seasonal differencing effectively removed the 24-hour cycle but left persistent 168-hour cyclicality in residual ACF functions. Introducing low-order Fourier terms into the regressor matrix eliminated weekly auto-dependencies without over-parameterizing the state space.

### 2. Volatility Filtering & Residual Diagnostics
Ljung-Box and ARCH-LM tests on squared residuals confirmed significant conditional heteroskedasticity in linear specifications. The asymmetric `eGARCH(1,1)-t` specification yielded the lowest AIC (**7.2485**). 

Post-estimation diagnostics confirmed that standardized residuals behaved as strict white noise, with an ARCH-LM test $p$-value of **0.8783** confirming the elimination of conditional volatility clusters.

### 3. Out-of-Sample Evaluation Strategy
* **Rolling Horizon Evaluation:** Expanding-window backtesting across the test partition re-estimating parameters every 24 hours.
* **24-Hour Horizon Forecast:** Out-of-sample forecast with dynamic 80% and 95% predictive density intervals evaluated against realized spot market prices.

---

## Performance & Benchmark Comparison

Point-forecast accuracy evaluated across the out-of-sample test horizon:

| Model Architecture | RMSE (EUR/MWh) | MAE (EUR/MWh) | Primary Advantage |
| :--- | :---: | :---: | :--- |
| **SARIMAX + eGARCH(1,1)-t** | **66.30** | **57.17** | Robust asymmetric volatility bounds & shock modeling |
| **Dynamic State Space (KFAS)** | **48.29** | **37.71** | Rapid parameter adaptation via dynamic latent states |

*Note: While the State Space model yielded lower point forecast error through dynamic state tracking, the SARIMAX-eGARCH model produced more reliable tail-risk predictive distributions.*

---

## Key Visualizations

### 1. SARIMAX Mean Equation Diagnostics
Residual diagnostics and ACF/PACF analysis verifying the elimination of autocorrelation and the effectiveness of harmonic Fourier regressors:
![SARIMAX Diagnostics](images/01_sarima_diagnostics.png)

### 2. eGARCH Volatility Modeling Diagnostics
Standardized residual diagnostics, QQ-plot, and ARCH-LM evaluation confirming the capture of conditional heteroskedasticity and fat tails:
![eGARCH Diagnostics](images/02_egarch_diagnostics.jpeg)

### 3. Pseudo Real-Time Backtesting
Expanding window evaluation tracking historical out-of-sample performance across recurring 24-hour trading horizons:
![Pseudo Real Time Forecasting](images/03_pseudo_real_time_forecasting.png)

### 4. SARIMAX-eGARCH 24-Hour Horizon Forecast
Out-of-sample day-ahead price trajectory accompanied by 80% and 95% dynamic predictive density intervals:
![Day-Ahead Forecast](images/04_sarimax_egarch_forecast.png)

### 5. Benchmark Performance: Econometric vs. State Space
Direct forecast trajectory comparison evaluating the hybrid SARIMAX-eGARCH model against the dynamic KFAS State Space baseline:
![Forecast Comparison](images/05_forecast_comparison_benchmark.jpeg)

---

## ⚙️ Repository Structure & Usage

```text
├── data/
│   ├── Finnish_Energy_Market_Weather_Data.csv     # Raw hourly prices and wind data
│   └── benchmark_state_space_forecast.csv          # State Space baseline estimates
├── docs/
│   └── finnish-electricity-price-forecasting.pdf
├── images/
│   ├── 01_sarima_diagnostics.png
│   ├── 02_egarch_diagnostics.jpeg
│   ├── 03_pseudo_real_time_forecasting.png
│   ├── 04_sarimax_egarch_forecast.png
│   └── 05_forecast_comparison_benchmark.jpeg
├── main.R                                          # Master R replication script
└── README.md

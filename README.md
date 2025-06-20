# Granger Causality & Stationarity via Fractional Differencing

This repository contains a small piece of code I wanted to preserve for future reference—and to share with anyone interested in the **Granger causality test** and the preprocessing required to apply it effectively.

**Note:** The code is not directly ready for production use; it's more of a conceptual prototype or "food for thought."

## Overview

The **Granger causality test** is an econometric hypothesis test used to determine whether one time-series variable provides statistically significant information to forecast another. It’s commonly applied to **multivariate time-series data** with a specified lag.

However, like many statistical tests in time-series analysis, the Granger test requires the input data to be **stationary**—that is, the series must have a constant mean, constant variance, and no seasonal structure.

## 🛠Methodology

This project proposes a preprocessing methodology to transform **non-stationary data into stationary data** using **fractional differencing**:

* Unlike traditional differencing (which may discard valuable information), **fractional derivatives** allow us to preserve more of the original signal from an **information-theoretic** perspective.
* The algorithm applies **incremental differencing with fractional orders**, followed by the **Augmented Dickey–Fuller (ADF)** test to check for stationarity.
* If the ADF test fails, the process continues iteratively until the series is deemed stationary.

## Notes

* If the series still fails to become stationary after second-order differencing, the issue may stem from **seasonality**, rather than a simple trend or variance shift. In such cases, the Granger causality test may not be applicable without additional seasonal adjustments.
* This idea is inspired by *Marcos López de Prado's* excellent book:
  **"Advances in Financial Machine Learning" (2018)** – see **Chapter 5: Fractionally Differentiated Features**. Highly recommended!

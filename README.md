This is a small piece of code that I wanted to keep here in case I need it later. I’m also sharing it for anyone interested in the Granger causality test and what it takes to use it.
The code is not directly usable; it’s just food for thought.

The Granger causality test is an econometric hypothesis test used to determine whether one variable helps forecast another in multivariate time‑series data with a particular lag.
As with many tests in the time‑series field, a prerequisite is that the data must be stationary, that is, it should have a constant mean, constant variance, and no seasonal component.

Here, I propose a methodology to transform non‑stationary data into stationary data by iteratively differencing it with fractional derivatives. I use fractional derivatives so that 
as little information as possible is lost from an information‑theoretic point of view. After applying these incremental differences with fractional coefficients, I run the augmented 
Dickey–Fuller (ADF) test to see whether the series has become stationary; if not, I keep looping.

Note: - If the data are still not stationary after second‑order differencing, the issue is probably a seasonal component rather than a non‑constant mean or variance. That’s the limit 
of this small project—in such a case, the Granger causality test cannot be used.
      - The idea comes from "Advances in Financial Machine Learning" (2018) by Marcos López de Prado - chapter 5 "Fractionally Differentiated Features" - check it out :)

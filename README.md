# Abstract

In survival analysis, right-censoring limits the observation of complete event times, making the validation of parametric assumptions challenging. This work presents a comparative assessment of three classical goodness-of-fit tests for right-censored data: Kolmogorov-Smirnov, Anderson-Darling, and Cramér-von Mises.

A Monte Carlo simulation study was conducted to evaluate how sample size, skewness, and censoring intensity affect test performance under administrative and random censoring mechanisms. The analysis considered Weibull, Lognormal, and Logistic distributions under composite null hypotheses, with procedures implemented using the `GofCens` R package and calibrated through bootstrap resampling.

The results show that administrative censoring strongly affects integrated distance-based statistics, making the Kolmogorov-Smirnov test the most reliable alternative in these scenarios. For random censoring, test performance depends on the interaction between data characteristics and censoring patterns. Additionally, this thesis introduces the $Diff_{t_{max}}$ diagnostic metric, which provides a practical indicator of tail distortion and supports the selection of appropriate goodness-of-fit procedures in censored survival datasets.

## Key Words
Survival Analysis, Right-Censoring, Goodness of fit, Monte Carlo Simulations, Test Selection

## Mathematics Subject Classification [MSC]
* **62N01:** Censored data models
* **62N03:** Testing
* **65C05:** Monte Carlo methods

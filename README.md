# HFMRD: Hedge Funds Misreported Returns Detector

This script represents a full-featured framework for detecting misreported returns in hedge funds. It be used in order to perform all the tests proposed by Bollen & Pool (2008):

* the `Primary Tests`: First Digits Analysis, Second Digits Analysis, First-Two Digits
* the `Advanced Tests`: Third Digits Analysis, Second Order Analysis, Summation Analysis
* the `Associated Tests`: Last-Two Digits Analysis, Number Duplication Analysis, Distortion Factor Model
* the `Mantissae Analysis`
* the `Zipf's Law Analysis`

For each significant digit analysis, the following conformity indicators are provided:

* **Goodness-of-Fit Measures (14):**
  * Anderson-Darling Discrete (Choulakian, 1994)
  * Chebyshev Distance (Leemis, 2000)
  * Cramer-von Mises Discrete (Choulakian, 1994)
  * Euclidean Distance (Cho & Gaines, 2007)
  * Freedman's U2 (Freedman, 1981)
  * Freeman-Tukey T2 (Freeman & Tukey, 1950)
  * Hotelling's Joint Digits (Hotelling, 1931)
  * Judge-Schechter Mean Deviation (Judge & Schechter, 2009)
  * Kolmogorov-Smirnov (Kolomonorgov, 1933)
  * Kuiper (Kuiper, 1960)
  * Likelihood Ratio (Neyman & Pearson, 1933)
  * Pearson's X2 (Pearson, 1900)
  * Watson's U2 Discrete (Choulakian, 1994)
* **Mean Absolute Deviation** (Nigrini et al., 2012)
* **Sum of Square Differences** (Kossovsky, 2014)
* **Z-Scores** (Nigrini et al., 2012)

## Requirements

The minimum Matlab version required is `R2014a`. In addition, the `Statistics and Machine Learning Toolbox` must be installed in order to properly execute the script.

## Dataset & Usage 

The framework doesn't require any specific dataset structure. Numeric data can be extracted from any source of produced following any methodology, but a minimum amount of 1000 elements (with at least 50 unique observations) is required in order to perform coherent analyses.

The `run.m` script provides an example of how this framework can be used, but all the functions located in the `Scripts` folder can be executed in standalone computation processes. It is recommended to validate and pre-process the dataset using the `benford_data` function. The `benford_analyse` functions can be used in order to perform a full automatic analysis of the dataset and plot the results. The `benford_random` function is an additional tool that produces random numbers whose digits follow the Benford's Law distribution.

## Screenshots

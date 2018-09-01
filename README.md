# HFMRD: Hedge Funds Misreported Returns Detector

This script represents a full-featured framework for detecting misreported returns in hedge funds through the following tests:

* `Low Correlation with Other Assets` (Bollen & Pool, 2008-2010):
  * IndexRSQ
  * MaxRSQ or SwitchRSQ with change-point regression (as per Andrews et al., 1996)
* `Serial Correlation` (Bollen & Pool, 2008-2010):
  * Unconditional
  * Conditional
* `Bias Ratio` (Abdulali, 2002)
* `December Spike` (Agarwal et al., 2011)
* `Discontinuity At Zero` (Bollen & Pool, 2008-2010)
* `Digits Conformity` (Bollen & Pool, 2008-2010):
  * Benford's Law Conformity of First Digits
  * Uniform Distribution Conformity of Last Digits
* `Data Quality` (Straumann, 2008):
  * Number of Negative Returns
  * Number of Zero Returns
  * Number of Unique Returns
  * Number of Pairs of Identical Returns
  * Maximum Length of Adjacent Identical Returns

## Requirements

The minimum Matlab version required is `R2014a`. In addition, the following products and toolboxes must be installed in order to properly execute the script:

* Financial Toolbox
* Image Processing Toolbox
* Mapping Toolbox
* Statistics and Machine Learning Toolbox

## Dataset & Usage 

The framework doesn't require any specific dataset structure. Numeric data can be extracted from any source of produced following any methodology, but a minimum amount of 1000 elements (with at least 50 unique observations) is required in order to perform coherent analyses.

The `run.m` script provides an example of how this framework can be used, but all the functions located in the `Scripts` folder can be executed in standalone computation processes. It is recommended to validate and pre-process the dataset using the `benford_data` function. The `benford_analyse` functions can be used in order to perform a full automatic analysis of the dataset and plot the results. The `benford_random` function is an additional tool that produces random numbers whose digits follow the Benford's Law distribution.

## Screenshots

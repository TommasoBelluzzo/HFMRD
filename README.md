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

## Usage

1. Create a properly structured database (see the paragraph below).
1. Edit the `run.m` script following your needs.
1. Execute the `run.m` script.

## Dataset

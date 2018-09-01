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

Every dataset must be structured like the default one included in any release of the framework (`Datasets/Example.xlsx`). The latter, based on the US financial sector, defines the following entities:

#### Benchmark (BM) & Risk-Free Rate (RF)

The benchmark is represented by the market proxy defined in Fama & French, 1993: the value-weighted returns of all the US CRSP firms listed on the AMEX, NASDAQ or NYSE that have a CRSP share code of 10 or 11 at the beginning of month t, good shares and price data at the beginning of t, and good return data for t. The 1M treasury bill rate is taken as the risk-free rate.

#### Hedge Funds (3):
* The Growth Fund of America - Class A (AGTHX)
* The Gateway Fund - Class A (GATEX)
* The Fairfield Sentry Fund of  Bernard Madoff (SENTRY)
		
#### Style Factors (18):
* **MRKEXC:** the excess return on the market, calculated as benchmark minus risk-free rate.
* **Fama & French Style Factors (8)** (taken from the [French Data Library](http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html))
  * American International Group Inc. (AIG)
  * The Allstate Corp. (ALL)
  * Berkshire Hathaway Inc. (BRK)
  * MetLife Inc. (MET)
  * Prudential Financial Inc. (PRU)
* **Group 2: Investment Banks (6)**
  * Bank of America Corp. (BAC)
  * Citigroup Inc. (C)
  * The Goldman Sachs Group Inc. (GS)
  * J.P. Morgan Chase & Co. (JPM)
  * Lehman Brothers Holdings Inc. (LEH)
  * Morgan Stanley (MS) 
* **Group 3: Commercial Banks (7)**
  * American Express Co. (AXP)
  * Bank of New York Mellon Corp. (BK)
  * Capital One Financial Corp. (COF)
  * PNC Financial Services Inc. (PNC)
  * State Street Corp. (STT)
  * US Bancorp (USB)
  * Wells Fargo & Co. (WFC)
* **Group 4: Government-sponsored Enterprises (2)**
  * American International Group Inc. (AIG)
  * The Allstate Corp. (ALL)
  
  
  










#### State Variables (6):

* **VIX:** the implied volatility index.
* **TBN3M:** the 3M treasury bill rate.
* **CRESPR:** the change in the credit spread (the BAA corporate bond rate minus the 10Y treasury bill rate).
* **LIQSPR:** the change in the liquidity spread (the 3M treasury bill rate minus the federal funds rate).
* **YIESPR:** the change in the yield spread (the 10Y treasury bond rate minus the 3M treasury bond rate).

For what concerns the financial time series:
* they must contain enough observations to run consistent calculations (a minimum of 252 observations for at least 5 firms is recommended);
* they must have been previously validated and preprocessed by discarding illiquid series with too many zeroes (unless necessary), detecting and removing outliers, removing rows with NaNs or filling the gaps with interpolation approach;
* returns must expressed on a logarithmic scale, in accordance with all the systemic risk indicators;
* market capitalizations must contain a supplementar observation at the beginning because a one-day lagged version is used in order to calculate weighted averages of probabilistic measures;
* total liabilities values must be rolled forward by at least three months (keeping a daily frequency) in order to simulate the difficulty of renegotiating debt in case of financial distress following the SRISK methodology;
* state variables and groups are optional, hence their respective sheets must be removed from the dataset if the related computations aren't necessary;
* state variables, if defined, must contain a supplementar observation at the beginning because a one-day lagged version is required in order to follow the CoVaR/ΔCoVaR methodology.

## Screenshots

![Probabilistic Measures](https://i.imgur.com/1Q1SQd2.png)

![Network Measures](https://i.imgur.com/NuSHgBO.png)

![Network Graph](https://i.imgur.com/fpEVHPf.png)


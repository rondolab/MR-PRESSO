# MR-PRESSO
MR-PRESSO (**Mendelian Randomization Pleiotropy RESidual Sum and Outlier**) is a framework that allows for the evaluation of pleiotropy in multi-instrument Mendelian Randomization utilizing genome-wide summary association statistics.

MR-PRESSO has three components, including:
1. detection of pleiotropy (*MR-PRESSO global test*)
2. correction of pleiotropy via outlier removal (*MR-PRESSO outlier test*)
3. testing of significant distortion in the causal estimates before and after outlier removal (*MR-PRESSO distortion test*).

### Reference

Widespread pleiotropy confounds causal relationships between complex traits and diseases inferred from Mendelian randomization. Marie Verbanck, Chia-Yen Chen, Benjamin Neale, Ron Do. bioRxiv 2017. DOI: 10.1101/157552.
<http://www.biorxiv.org/content/early/2017/06/30/157552>

### 1. Install and load MR-PRESSO
To install the latest development builds directly from GitHub, run this instead:
```r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
```
Load MR-PRESSO 
```r
library(MRPRESSO)
```

### 2. Example
```r
# Load a simulated toy dataset
data(SummaryStats)

# Run MR-PRESSO global framework
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)
```

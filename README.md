# MR-PRESSO
MR-PRESSO (**Mendelian Randomization Pleiotropy RESidual Sum and Outlier**) is a method that allows for the evaluation of horizontal pleiotropy in multi-instrument Mendelian Randomization utilizing genome-wide summary association statistics.

MR-PRESSO has three components, including:
1. detection of horizontal pleiotropy (*MR-PRESSO global test*)
2. correction of horizontal pleiotropy via outlier removal (*MR-PRESSO outlier test*)
3. testing of significant distortion in the causal estimates before and after outlier removal (*MR-PRESSO distortion test*).

### Reference

Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. Marie Verbanck, Chia-Yen Chen, Benjamin Neale, Ron Do. Nature Genetics 2018. DOI: 10.1038/s41588-018-0099-7.
<https://www.nature.com/articles/s41588-018-0099-7>

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

# Run MR-PRESSO global method
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)
```

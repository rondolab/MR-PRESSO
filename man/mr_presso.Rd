\name{mr_presso}
\alias{mr_presso}
\title{
    a function to perform the Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO) test
}
\description{
    MR-PRESSO (Mendelian Randomization Pleiotropy RESidual Sum and Outlier) is a unified framework that allows for the evaluation of pleiotropy in a standard MR model. The method extends on previous approaches that utilize the general model of multi-instrument MR on summary statistics. MR-PRESSO has three components, including: 1) detection of pleiotropy (MR-PRESSO global test); 2) correction of pleiotropy via outlier removal (MR-PRESSO outlier test); and 3) testing of significant distortion differences in the causal estimates before and after MR-PRESSO correction (MR-PRESSO bias assessment test).
}
\usage{
    mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, OUTLIERtest = FALSE, BIAStest = FALSE, SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL)
}
\arguments{
  \item{BetaOutcome}{
    character, name of the outcome variable
}
  \item{BetaExposure}{
    vector of characters, name(s) of the exposure(s)
}
  \item{SdOutcome}{
    character, name of the standard deviation of the outcome variable
}
  \item{SdExposure}{
    vector of characters, name(s) of the standard deviation of the exposure(s)
}
  \item{data}{
    dataframe of summary statistics containing the outcome and exposure variables 
}
  \item{OUTLIERtest}{
    boolean, if TRUE the MR-PRESSO outlier test will be performed, default is FALSE
}
  \item{BIAStest}{
    boolean, if TRUE the MR-PRESSO bias test on the causal estimate will be performed, default is FALSE
}
  \item{SignifThreshold}{
    float, significance threshold to use between 0 and 1, default is 0.05
}
  \item{NbDistribution}{
    integer, number of elements to simulate to form the null distribution to compute empirical P-values. This is directly impacting the precision of the empirical P-values, it is 1/NbDistribution for the global test and nrow(data)/NbDistribution
}
  \item{seed}{
    a single value, interpreted as an integer to use in set.seed()
}

}

% \details{}

\value{
    \item{Global Test}{Results of the MR-PRESSO global test. List of the observed residual sum of squares 'RSSobs' and empirical P-value 'Pvalue'}
    \item{Outlier Test }{Results of the MR-PRESSO outlier test. Table of observed residual sum of squares 'RSSobs' and 'Pvalue' per SNV}
    \item{Bias Test}{Results of the MR-PRESSO bias assessment test. List of the 'Outliers Indices' identified as outliers and excluded to calculate the bias as the 'Distortion Coefficient' and its 'Pvalue'.}
}
\references{
    Widespread pleiotropy confounds causal relationships between complex traits and diseases inferred from Mendelian randomization. Marie Verbanck*, Chia-Yen Chen*, Benjamin Neale, Ron Do.
}
\author{
    Marie Verbanck
}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{}

\examples{

data(SummaryStats)
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, BIAStest = TRUE, data = SummaryStats, NbDistribution = 100,  SignifThreshold = 0.05)
}
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

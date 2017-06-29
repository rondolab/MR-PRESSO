\name{MRPRESSO-package}
\alias{MRPRESSO-package}
\alias{MRPRESSO}
\docType{package}
\title{
\packageTitle{MRPRESSO}
}
\description{
MR-PRESSO (Mendelian Randomization Pleiotropy RESidual Sum and Outlier) is a unified framework that allows for the evaluation of pleiotropy in a standard MR model. The method extends on previous approaches that utilize the general model of multi-instrument MR on summary statistics. MR-PRESSO has three components, including: 1) detection of pleiotropy (MR-PRESSO global test); 2) correction of pleiotropy via outlier removal (MR-PRESSO outlier test); and 3) testing of significant distortion differences in the causal estimates before and after MR-PRESSO correction (MR-PRESSO bias assessment test).
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MRPRESSO}
\packageIndices{MRPRESSO}
}
\author{
\packageAuthor{MRPRESSO}

Maintainer: \packageMaintainer{MRPRESSO}
}
\references{
    Widespread pleiotropy confounds causal relationships between complex traits and diseases inferred from Mendelian randomization. Marie Verbanck*, Chia-Yen Chen*, Benjamin Neale, Ron Do.
}
% \keyword{ package }
% \seealso{}
\examples{
data(SummaryStats)
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, BIAStest = TRUE, data = SummaryStats, NbDistribution = 100,  SignifThreshold = 0.05)
}

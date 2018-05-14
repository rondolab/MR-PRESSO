\name{SummaryStats}
\alias{SummaryStats}
\docType{data}
\title{
Simulated toy dataset
}
\description{
'SummaryStats' is a simulated toy dataset of summary statistics for Y the outcome and E1 and E2 two exposures. It is composed of 45 single ncleotide variants (SNVs) associated with E1 and 5 additional variants associated with E1 and E2 which are therefore submitted to pleiotropy. Test the function mr_presso() on this toy dataset.
}
\usage{data("SummaryStats")}
\format{
  A data frame with 50 observations on the following 9 variables.
  \describe{
    \item{\code{E1_effect}}{a numeric vector}
    \item{\code{E1_se}}{a numeric vector}
    \item{\code{E1_pval}}{a numeric vector}
    \item{\code{E2_effect}}{a numeric vector}
    \item{\code{E2_se}}{a numeric vector}
    \item{\code{E2_pval}}{a numeric vector}
    \item{\code{Y_effect}}{a numeric vector}
    \item{\code{Y_se}}{a numeric vector}
    \item{\code{Y_pval}}{a numeric vector}
  }
}
% \details{}
% \source{}

\references{
    Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. Marie Verbanck, Chia-Yen Chen, Benjamin Neale, Ron Do. Nature Genetics 2018. DOI: 10.1038/s41588-018-0099-7. https://www.nature.com/articles/s41588-018-0099-7
}
\examples{
data(SummaryStats)
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)
}
\keyword{datasets}

\name{BiomarkeR-package}
\alias{BiomarkeR-package}
\alias{BiomarkeR}
\docType{package}
\title{
Paired and Unpaired Biomarker Identifier including a method to infer networks.
}
\description{
This package allows to rank features and infer networks based on the Paired (pBI) and Unpaired Biomarker Identifier (uBI). 
}
\details{
\tabular{ll}{
Package: \tab BiomarkeR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-26\cr
License: \tab LGPL \cr
LazyLoad: \tab yes\cr
}
}
\author{
Michael Netzer <michael.netzer@umit.at>
}
\references{
Baumgartner, C., Lewis, G. D., Netzer, M., Pfeifer, B., and Gerszten, R. E. (2010). A new data mining approach for profiling and categorizing kinetic patterns of
metabolic biomarkers after myocardial injury. Bioinformatics, 26(14), 1745-751
}

\keyword{biomarker discovery}
\seealso{
}
\examples{
data(BIdata)
pbi.scores <- pBI(dataset=pbi.data, classlabels=pbi.class, referenceclasslabel="rest", ids=pbi.id, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 5)
ubi.scores <- uBI(dataset=ubi.data, classlabels=ubi.class, referenceclasslabel="control", useMedian = TRUE, lambda = 100, plotScores = FALSE, numTopRankedToPlot = 5)
}

\name{uBI}
\alias{uBI}
\title{
Unpaired Biomarker Identifier (uBI)
}
\description{
Calculates a score for every attribute representing its discriminatory ability. The uBI method adresses independent samples.
}
\usage{
uBI(dataset, classlabels, referenceclasslabel, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
	Dataset where rows are attributes and columns are samples.
}
  \item{classlabels}{
	A vector of classlabels. Number of different classlabels must be 2.
}
  \item{referenceclasslabel}{
	Classlabel of the reference class.
}
  \item{useMedian}{
	Whether to use median (default) or mean as location parameter.
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{plotScores}{
	If true, a barplot with top ranked attributes is created.
}
  \item{numTopRankedToPlot}{
	Number of top ranked attributes to be plotted (if plotScores=TRUE).
}
}
\details{
Unpaired Biomarker Identifier (uBI) calculates a score for every attribute representing its discriminatory ability adressing independent samples. 
The score is based on the product of sensitivity and specificity (TP*), the relative change compared to the reference group (delta) and the ratio of the 
coeffiecients of variation in the referenceversus the comparison group (CVref/CV). 
}
\value{
A vector, containing the pBI scores. 
}
\references{
Baumgartner, C., Lewis, G. D., Netzer, M., Pfeifer, B., and Gerszten, R. E. (2010). A new data mining approach for profiling and categorizing kinetic patterns of
metabolic biomarkers after myocardial injury. Bioinformatics, 26(14), 1745-1751
}
\author{
Michael Netzer, michael.netzer@umit.at
}
\note{
}

\seealso{\code{\link{pBI}}, \code{\link{pBIGraph}}, \code{\link{uBIGraph}}}

\examples{
data(BIdata)
ubi.scores <- uBI(dataset=ubi.data, classlabels=ubi.class, referenceclasslabel="control", useMedian = TRUE, lambda = 100, plotScores = FALSE, numTopRankedToPlot = 5)
}
\keyword{ univar }

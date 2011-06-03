\name{pBI}
\alias{pBI}

\title{Paired Biomarker Identifier (pBI)}
\description{Calculates a score for every attribute representing its discriminatory ability. The pBI method adresses dependent samples.}
\usage{
pBI(dataset, classlabels, referenceclasslabel, ids, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 10, bars.cols = c("red","blue"))
}
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
  \item{ids}{
	A vector, containing the ids for the samples. 
}
  \item{useMedian}{
	Whether to use median (default) or mean as location parameter.
}
  \item{lambda}{
	Scaling factor for the scores. 
}
  \item{plotScores}{
	If true, a barplot with top ranked attributes is created.
}
  \item{numTopRankedToPlot}{
	Number of top ranked attributes to be plotted (if plotScores=TRUE).
}
  \item{bars.cols}{
	Vector with two colors (first color for bars >= 0, second color for bars < 0).
}
}
\details{
Paired Biomarker Identifier (pBI) calculates a score for every attribute representing its discriminatory ability adressing dependent samples. 
The score is based on the discriminance measure (DA), the median percent change (delta) and the coeffiecient of variation (CV). 
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

\seealso{\code{\link{uBI}}, \code{\link{pBIGraph}}, \code{\link{uBIGraph}}}

\examples{
data(BIdata)
pbi.scores <- pBI(dataset=pbi.data, classlabels=pbi.class, referenceclasslabel="rest", ids=pbi.id, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ univar }

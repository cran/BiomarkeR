\name{uBIGraph}
\alias{uBIGraph}

\title{
Unpaired Biomarker Identifier (uBI) graph.
}
\description{
Infers a graph based on the uBI scores (adresses independent samples). 
}
\usage{
uBIGraph(dataset, classlabels, referenceclasslabel, useMedian = TRUE, lambda = 100, threshold = "q90", plotGraph = FALSE, edge.file = NULL)
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
	Scaling factor for the scores. 
}
  \item{threshold}{
	Threshold for defining an edge in the graph. A edge between vertex X (metabolite X) and vertex Y (metabolite Y) is drawn if the
	score of the ratio is greater than this threshold. Threshold is selected as the 90\% percentile of scores by default. 
}
  \item{plotGraph}{
	If true, the graph is plotted. 
}
  \item{edge.file}{
	An edge file is saved into the given path (if edge.file!=NULL).
}
}
\details{
	Infers a graph based on the pBI scores. Vertices in such networks represent attributes (metabolites).
	An edge between vertex X (metabolite X) and vertex Y (metabolite Y) is drawn if the
	score of the ratio is greater than a user defined threshold. 
}
\value{
A graph as a igraph object.
}
\references{
}
\author{
Michael Netzer, michael.netzer@umit.at
}
\note{
}

\seealso{\code{\link{pBIGraph}}, \code{\link{pBI}}, \code{\link{uBI}}}

\examples{
data(BIdata)
g <- uBIGraph(dataset=ubi.data, classlabels=ubi.class, referenceclasslabel="control", useMedian = TRUE, lambda = 100, threshold = "q90", plotGraph = FALSE, edge.file = NULL)
}
\keyword{ multivariate }

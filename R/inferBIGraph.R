pBIGraph <- function(dataset, classlabels, referenceclasslabel, ids, useMedian = TRUE, lambda = 100, threshold = "q90", plotGraph = FALSE, edge.file = NULL)
{
	if(any(dataset==0))
		stop("No zero values are allowed!")
	ratios  <- .calculateRatios(dataset)
	biscores <- pBI(ratios, classlabels, referenceclasslabel, ids, useMedian = useMedian, lambda = lambda, plotScores = FALSE)
	g <- .createGraph(biscores, nodes= rownames(dataset), threshold, edge.file, plotGraph)
	return(g)
}

uBIGraph <- function(dataset, classlabels, referenceclasslabel, useMedian = TRUE, lambda = 100, threshold = "q90", plotGraph = FALSE, edge.file = NULL)
{
	if(any(dataset==0))
		stop("No zero values are allowed!")
	ratios  <- .calculateRatios(dataset)
	biscores <- uBI(ratios, classlabels, referenceclasslabel, useMedian = useMedian, lambda = lambda, plotScores = FALSE)
	g <- .createGraph(biscores, nodes= rownames(dataset), threshold, edge.file, plotGraph)
	return(g)
}

.createGraph <- function(biscores, nodes, threshold, edge.file, plotGraph)
{
	significant <- names(.getSignificantValues(abs(biscores), threshold))
	g <- .getGraph(nodes, significant)
	if(!is.null(edge.file))
		{
			edgelist <- .getEdgelist(nodes, significant)
			write.table(edgelist,file=edge.file, row.names = FALSE, col.names = FALSE, quote=FALSE)
		}
	if(plotGraph)
		plot(g)
	return(g)
}

.calculateRatios <- function(dataset){
	temp <- NULL
	tnames <- NULL
	for(i in 1:(nrow(dataset)-1)){
		for(j in (i+1):nrow(dataset))
		{
			temp <- rbind(temp, na.keep(dataset[i,]/dataset[j,]))
			tnames <- c(tnames, paste(rownames(dataset)[i],rownames(dataset)[j], collapse="", sep="/"))
		}
	}
	temp <- data.frame(temp)  
	rownames(temp) <- tnames
	return(as.matrix(temp))
}

.getSignificantValues <- function(scores, threshold)
{
	if(substr(threshold,1,1)=="q")
		{
			quant <- substr(threshold, 2, nchar(threshold))
			if(is.na(as.numeric(quant)))
				stop("Quantile must be numeric!")
			threshold <- quantile(scores, as.numeric(quant)/100)
		}
	else if(!is.na(as.numeric(threshold)))
			threshold <- as.numeric(threshold)
	else
		stop("Error: Threshold must be either a quantile or a numerical value!")
	significant <- scores[scores > threshold]
	return(significant)
}

.getEdgelist <- function(nodes, edges, edge.splitchar = "/")
{
	edgelist <- matrix(data=NA, ncol=2, nrow=length(edges))
	for(i in 1:length(edges)){
	   temp <- edges[i]
	   temp <- unlist(strsplit(temp, edge.splitchar))
	   edgelist[i,] <- temp
	}
	return(edgelist)
}

.getGraph <- function(nodes, edges, edge.splitchar = "/")
{
	if(length(nodes)==0)
		stop("Number of nodes is zero!")
	g <- NULL
	g <- new("graphNEL", nodes = nodes)
	for(i in 1:length(edges)){
	   temp <- edges[i]
	   temp <- unlist(strsplit(temp, edge.splitchar))
	   g <- addEdge(temp[1], temp[2], g)
	}
	return(g)
}
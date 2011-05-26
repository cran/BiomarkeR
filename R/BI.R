pBI <- function(dataset, classlabels, referenceclasslabel, ids, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 10)
{
	if(!is.matrix(dataset))
		stop("Dataset must be a matrix!")
	if(!is.vector(classlabels))
		stop("Classlabel must be a vector!")
	if(!is.vector(ids))
		stop("Ids must be a vector!")
	if(length(unique(classlabels))!=2)
		stop("Number of class labels must be two!")
	if(!referenceclasslabel%in%classlabels)
		stop("Reference classlabel must be a valid classlabel!")
	if(length(ids)!=ncol(dataset))
		stop("Number of ids must be equal to number of samples!")
	if(sum(is.na(dataset))>0)
		stop("No missing values are allowed!")
	if(any(dataset==0))
		stop("No zero values are allowed!")
	dataset.reference <-  dataset[,classlabels==referenceclasslabel]
	dataset.comp <- dataset[,classlabels!=referenceclasslabel]
	dataset.reference.ids <- ids[classlabels==referenceclasslabel]
	dataset.comp.ids <- ids[classlabels!=referenceclasslabel]
	dataset.comp <- dataset.comp[, match(dataset.reference.ids, dataset.comp.ids)]
	if(sum(dim(dataset.comp)!=dim(dataset.reference))> 0)
		stop("Dataset to compare and reference dataset must have same dimension!")
	dataset <- .calcPercentageChange(dataset.comp, dataset.reference)

	 pBIScores <- apply(dataset, 1, function(x){
		da <- .getDA(x, mu=0)
		locPar <- ifelse(useMedian, median(x), mean(x))
		delta <- ((abs(locPar)/100) + 1 ) * sign(locPar)
		if(mean(x) == 0) 
			pBI <- 0
		else {
				cv <- sd(x)/mean(x)
				if(abs(cv) > 1)
					cv <- 1
				pBI <- lambda * da *  sqrt(abs(delta/cv))* sign(delta)
			}
		return(pBI)
	 })
	 names(pBIScores) <- rownames(dataset)
	if(plotScores)
		.plotScores(pBIScores, numTopRankedToPlot, method="pBI")
	return(pBIScores)
}

uBI <- function(dataset, classlabels, referenceclasslabel, useMedian = TRUE, lambda = 100, plotScores = TRUE, numTopRankedToPlot = 10)
{
	if(!is.matrix(dataset))
		stop("Dataset must be a matrix!")
	if(!is.vector(classlabels))
		stop("Classlabel must be a vector!")
	if(length(unique(classlabels))!=2)
		stop("Number of class labels must be two!")
	if(!referenceclasslabel%in%classlabels)
		stop("Reference classlabel must be a valid classlabel!")
	if(sum(is.na(dataset))> 0)
		stop("No missing values are allowed!")
	if(any(dataset==0))
		stop("No zero values are allowed!")
	uBIScores <- apply(dataset, 1, function(x){
		tp2 <- .getTP2(x, classlabels)
		x.ref <- x[classlabels==referenceclasslabel]
		x.comp <- x[classlabels!=referenceclasslabel]
		if(useMedian){
				locPar.ref <- median(x.ref)
				locPar.comp <- median(x.comp)
			}
		else {
				locPar.ref <- mean(x.ref)
				locPar.comp <- mean(x.comp)
			}
		if(locPar.ref != 0 & mean(x.comp)!=0 & mean(x.ref)!=0){ 
				fDiff = locPar.comp/locPar.ref
				delta <- ifelse(fDiff >= 1, fDiff,  -(1/fDiff))
				cv.comp <- sd(x.comp)/mean(x.comp)
				cv.ref <- sd(x.ref)/mean(x.ref)
				ratio = cv.ref / cv.comp
				if (ratio < 1)
						ratio <- 1
				uBI <- lambda * tp2 * sqrt(abs(ratio * delta))* sign(delta)
			}
		else 
			uBI <- Inf
		
		return(uBI)
	 })
	 names(uBIScores) <- rownames(dataset)
	if(plotScores)
		.plotScores(uBIScores, numTopRankedToPlot, method="uBI")
	return(uBIScores)
}

.plotScores <- function(scores, numTopRankedToPlot, method)
{
		numTopRanked <- min(numTopRankedToPlot, length(scores))
		scoresRanked <- scores[match(sort(abs(scores), decreasing = TRUE), abs(scores))]
		scoresTopRanked <- scoresRanked[1:numTopRanked]
		cols <- ifelse(scoresTopRanked >= 0, "red", "blue")
		zoom <- 1.5
		barplot(scoresTopRanked, main =paste(method, " scores (top ", numTopRanked, "ranked attributes)", sep=" "), ylab=method, col=cols, cex.axis = zoom, cex.names=zoom, cex.lab=zoom, cex.main=zoom)  
}

.getTP2<- function (values, classlabels) {
		 tmp <- as.data.frame(cbind(values, classlabels))
		 tmp$classlabels <- as.factor( tmp$classlabels)
		 model <- glm(classlabels~values, data=tmp, family=binomial(link = "logit"))
		confMat <- table(fitted(model)> 0.5, tmp$classlabels)
		tp2 <- 1
		for(i in 1:length(unique(classlabels)))
		{
			tp <- (diag(confMat)/apply(confMat, 2, sum))[i]
			if(tp < 0.5)
				tp <- 0
			tp2 <- tp2 * tp
		}
		tp2 <- ((tp2-0.25)*4/3)
		return(tp2)		
}

.getDA <- function (values, mu) {
		 n <- length(values)
		halfN <- n / 2
		numNeutral <- sum(values == mu)
		numUpper <- sum(values > mu)
		numLower <- sum(values < mu)
		numMax <- max(numUpper, numLower)
		numAcceptance = numMax + numNeutral
		return ((numAcceptance / halfN) - 1)
}

.calcPercentageChange <- function(dataset, reference)
{
	if(sum(dim(dataset)!=dim(reference))> 0)
		stop("Dataset to compare and reference dataset must have same dimension!")
	percMat <- matrix(data=NA, nrow = nrow(dataset), ncol = ncol(dataset), dimnames = dimnames(dataset))
	for(i in 1:ncol(dataset))
	{	
		ratios <- dataset[,i]/reference[,i]
		percentages <- vector()
		for(j in 1:length(ratios))
		{
			tmp <- ratios[j]
			if (tmp >= 1)
				percentages[j] <- ((tmp-1)*100)
			else 
				percentages[j] <- ((1/tmp)-1)   *-100
		}
		percMat[,i] <- percentages
	}
	return(percMat)
}
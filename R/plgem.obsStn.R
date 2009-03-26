"plgem.obsStn" <-
function(data, plgemFit, covariate=1, baselineCondition=1, verbose=FALSE) {

	#some checks...
	.checkExpressionSet(data)
	if(class(plgemFit)!="list") stop("Object plgemFit in function plgem.obsStn is not of class list")

  covariate <- .checkCovariate(covariate, pData(data))
	condition.names <- as.character(pData(data)[, covariate])

  baselineCondition <- .checkCondition(baselineCondition, "baselineCondition", covariate, pData(data))

	if(class(verbose)!="logical") stop("Argument verbose in function plgem.obsStn is not of class logical")
	
	if(verbose) cat("calculating observed PLGEM-STN statistics:")

	# preparing...
	condition.names <- as.character(pData(data)[, covariate])
	condition.name <- unique(condition.names)
	condition.number <- length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object data for function plgem.obsStn")
	if(verbose) cat("found",(condition.number-1),"condition(s) to compare to the baseline.\n")
	dataMatrix <- exprs(data)
	
	# replacing zero and negative values with minimum positive value
	dataMatrix <- replace(dataMatrix, dataMatrix<=0, min(dataMatrix[dataMatrix > 0]))
	if(verbose) cat("working on baseline", baselineCondition, "...\n")
	baseline.col <- which(condition.names == baselineCondition)
	if (verbose) cat(colnames(dataMatrix)[baseline.col], "\n")
	observedStn <- array(, dim=c(nrow(dataMatrix), condition.number - 1))
	rownames(observedStn) <- featureNames(data)
	colnames(observedStn) <- condition.name[-which(condition.name == baselineCondition)]

	#calculating mean and modeled spread for the baseline condition
	mean.left <- rowMeans(dataMatrix[, baseline.col], na.rm=TRUE)
	spread.left <- .plgemSpread(mean.left, plgemFit$SLOPE, plgemFit$INTERCEPT)

	#calculating mean and modeled spread for the remaining condition(s)
	col.counter <- 0
	for (i in (1:condition.number)[-which(condition.name == baselineCondition)]) {
		col.counter <- col.counter + 1
		if(verbose) cat("working on condition", condition.name[i], "...\n")
		condition.col<-which(condition.names == condition.name[i])
		if(verbose) cat(colnames(dataMatrix)[condition.col], "\n")
		if(length(condition.col) == 1) {
      mean.right<-dataMatrix[, condition.col]
    }	else {
		  mean.right<-rowMeans(dataMatrix[, condition.col], na.rm=TRUE)
		}
		spread.right <- .plgemSpread(mean.right, plgemFit$SLOPE, plgemFit$INTERCEPT)
		#computation of PLGEM-STN statistics
		observedStn[, col.counter] <- .stn(mean.left, mean.right, spread.left, spread.right)
	}

	if(verbose) cat("done with calculating PLGEM-STN statistics.\n\n")
	gc()
	return(list("fit"=plgemFit, "PLGEM.STN"=observedStn))
}

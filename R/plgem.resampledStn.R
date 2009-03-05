"plgem.resampledStn" <-
function(data, plgemFit, covariate=1, baselineCondition=1, iterations="automatic", verbose=FALSE) {

	#some checks...
	.checkExpressionSet(data)

  covariate <- .checkCovariate(covariate, pData(data))
	condition.names <- as.character(pData(data)[, covariate])

  baselineCondition <- .checkCondition(baselineCondition, "baselineCondition", covariate, pData(data))

	if(class(plgemFit)!="list") stop("Object plgemFit in function plgem.resampledStn is not of class list")
	if(iterations!="automatic" && class(iterations)!="numeric" && class(iterations)!="integer") stop("Argument iterations in function plgem.resampledStn is neither of class numeric (or integer) nor equal to 'automatic'")
	if(class(verbose)!="logical") stop("Argument verbose in function plgem.resampledStn is not of class logical")

	if(verbose) cat("calculating resampled PLGEM-STN statistics:")

	# preparing...
	condition.name<-unique(condition.names)
	condition.number<-length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object data for function plgem.resampledStn")
	if (verbose) cat("found",(condition.number-1),"condition(s) to compare to the baseline.\n")
	dataMatrix<-exprs(data)
	
	# replacing zero and negative values with minimum positive value
	dataMatrix<-replace(dataMatrix,dataMatrix<=0,min(dataMatrix[dataMatrix>0]))
	rowNumber<-nrow(dataMatrix)
	baseline.col<-which(condition.names == baselineCondition)
	if (verbose) cat("baseline samples:\n")
	if (verbose) cat(colnames(dataMatrix)[baseline.col],"\n")
	fit.col <- which(condition.names == plgemFit$FIT.CONDITION)
	if (verbose) cat("resampling on samples:\n")
	if (verbose) cat(colnames(dataMatrix)[fit.col],"\n")

	#determining the number of replicates of each condition
	repl.number <- table(condition.names)
  repl.cases <- unique(
    repl.number[-which(names(repl.number) == baselineCondition)]
  )
  
	#determination of number of iterations
	if(iterations == "automatic") {
		a<-length(fit.col)
		b<-length(baseline.col)
		c<-max(repl.number[-1])
		iterations<-a^(b+c)
		iterations<-min(iterations, 500)
	}
	if(verbose) cat("Using", iterations, "iterations...\n")

	#computing resampled STN statistics for each case of number of replicates
	resampledStn<-array(,dim=c(rowNumber*iterations,length(repl.cases)))
	colnames(resampledStn)<-as.character(repl.cases)
	for (i in 1:length(repl.cases)) {
		if(verbose) cat("working on cases with", repl.cases[i], "replicates...\n")
		if(verbose) cat("     Iterations: ")
		for (j in 1:iterations){
			if(verbose) {if (j/20 == trunc(j/20)) {cat(j," ")} }
			#sampling column indices
			left.col<-sample(fit.col,length(baseline.col),replace=TRUE)
			right.col<-sample(fit.col,repl.cases[i],replace=TRUE)

			#calculating mean and modeled spread for the first artificial condition
			mean.left<-rowMeans(dataMatrix[,left.col],na.rm=TRUE)
			spread.left<-.plgemSpread(mean.left,plgemFit$SLOPE,plgemFit$INTERCEPT)	   

			#calculating mean and modeled spread for the second artificial condition
			if(length(right.col)==1) {
        mean.right <- dataMatrix[,right.col]
      } else {
			  mean.right <- rowMeans(dataMatrix[, right.col], na.rm=TRUE)
			}
			spread.right <- .plgemSpread(mean.right, plgemFit$SLOPE, plgemFit$INTERCEPT)

			#computation of resampled PLGEM-STN statistics
			rowIndex<-(rowNumber*(j-1)+1):(rowNumber*j)
			resampledStn[rowIndex,i]<-.stn(mean.left,mean.right,spread.left,spread.right)
		}
		if(verbose) cat("\n")
	}
	if(verbose) cat("done with calculating resampled PLGEM-STN statistics.\n\n")
	gc()
	return(list(RESAMPLED.STN=resampledStn, REPL.NUMBER=repl.number))
}

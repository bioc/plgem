"plgem.obsStn" <-
function(data,plgemFit,baseline.condition=1,verbose=FALSE) {
	library(Biobase)

	#some checks...
	if(class(data)!="exprSet") stop("Object data in function plgem.obsStn is not of class exprSet")
	if(!("conditionName" %in% colnames(pData(data)))) stop("the covariate conditionName is not defined in the data exprSet")
	if(class(plgemFit)!="list") stop("Object plgemFit in function plgem.obsStn is not of class list")
	if(class(baseline.condition)!="numeric" && class(baseline.condition)!="integer") stop("Argument baseline.condition in function plgem.obsStn is not of class numeric or integer")
	if(class(verbose)!="logical") stop("Argument verbose in function plgem.obsStn is not of class logical")
	
	if(verbose) cat("calculating observed PLGEM-STN statistics:")

	#internal functions
	stn<-function(location1,location2,spread1,spread2){
		(location2-location1)/(spread1+spread2)
	}

	plgem.spread<-function(location,slope,intercept) {
		(location^(slope))*exp(intercept)
	}
	
	#preparing...
	condition.name<-unique(as.character(data$conditionName))
	condition.number<-length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object data for function plgem.obsStn")
	if(verbose) cat("found",(condition.number-1),"condition(s) to compare to the baseline.\n")
	dataMatrix<-exprs(data)
	#replacing zero and negative values with minimum positive value
	dataMatrix<-replace(dataMatrix,dataMatrix<=0,min(dataMatrix[dataMatrix>0]))
	if(verbose) cat("working on baseline",condition.name[baseline.condition],"...\n")
	baseline.col<-which(data$conditionName==condition.name[baseline.condition])
	if (verbose) cat(colnames(dataMatrix)[baseline.col],"\n")
	obervedStn<-array(,dim=c(nrow(dataMatrix),condition.number-1))
	rownames(obervedStn)<-geneNames(data)
	colnames(obervedStn)<-condition.name[-baseline.condition]

	#calculating mean and modeled spread for the baseline condition
	mean.left<-rowMeans(dataMatrix[,baseline.col],na.rm=TRUE)
	spread.left<-plgem.spread(mean.left,plgemFit$SLOPE,plgemFit$INTERCEPT)

	#calculating mean and modeled spread for the remaining condition(s)
	col.counter<-0
	for (i in (1:condition.number)[-baseline.condition]) {
		col.counter<-col.counter+1
		if(verbose) cat("working on condition",condition.name[i],"...\n")
		condition.col<-which(data$conditionName==condition.name[i])
		if(verbose) cat(colnames(dataMatrix)[condition.col],"\n")
		if(length(condition.col)==1) {mean.right<-dataMatrix[,condition.col]}
		else {
		    mean.right<-rowMeans(dataMatrix[,condition.col],na.rm=TRUE)
		}
		spread.right<-plgem.spread(mean.right,plgemFit$SLOPE,plgemFit$INTERCEPT)
		#computation of PLGEM-STN statistics
		obervedStn[,col.counter]<-stn(mean.left,mean.right,spread.left,spread.right)
	}

	if(verbose) cat("done with calculating PLGEM-STN statistics.\n\n")
	gc()
	return(obervedStn)
}

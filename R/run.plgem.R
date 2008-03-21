"run.plgem" <-
function(esdata, signLev=0.001, rank=100, covariateNumb=1, baselineCondition=1, Iterations="automatic", fitting.eval=TRUE, plotFile=FALSE, writeFiles = FALSE, Verbose=FALSE) {

	#some checks
	if(class(esdata)!="ExpressionSet") stop("Argument 'esdata' is not of class 'ExpressionSet'.")
	if(class(signLev)!="numeric" && class(signLev)!="integer") stop("Argument 'signLev' is not of class 'numeric' or 'integer'.")
	if(class(rank)!="numeric" && class(rank)!="integer") stop("Argument 'rank' is not of class 'numeric' or 'integer'.")
	if(class(covariateNumb)!="numeric" && class(covariateNumb)!="integer") stop("Argument 'covariateNumb' is not of class 'numeric' or 'integer'.")
  if(as.integer(covariateNumb) > ncol(pData(esdata))) stop("Argument 'covariateNumb' is greater than the number of covariates in 'esdata'.")
	if(class(baselineCondition)!="numeric" && class(baselineCondition)!="integer") stop("Argument 'baselineCondition' is not of class 'numeric' or 'integer'.")
	if(Iterations!="automatic" && class(Iterations)!="numeric" && class(Iterations)!="integer") stop("Argument 'Iterations' is neither of class 'numeric' (or 'integer') nor equal to 'automatic'.")
	if(class(fitting.eval)!="logical") stop("Argument 'fitting.eval' is not of class 'logical'.")
	if(class(plotFile)!="logical") stop("Argument 'plotFile' is not of class 'logical'.")
	if(class(Verbose)!="logical") stop("Argument 'Verbose' is not of class 'logical'.")

  covariateNumb <- as.integer(covariateNumb)
  condition.names<-as.character(pData(esdata)[,covariateNumb])
	condition.name<-unique(condition.names)
	condition.number<-length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object 'esdata' for function 'run.plgem'.")

  baselineCondition <- as.integer(baselineCondition)
	if (baselineCondition > condition.number) stop("Argument 'baselineCondition' is greater than the number of conditions in 'esdata'.")

	#determining the number of replicates of each condition
	repl.number<-array(,dim=length(condition.name))
	for (i in 1:length(condition.name)) {
	    repl.number[i]<-length(which(condition.names==condition.name[i]))
	}
	if(max(repl.number)==1) stop ("PLGEM can not be fitted without replicates.")
	names(repl.number)<-condition.name

	#determination of the best condition on which to fit the model
	if(length(which(repl.number==max(repl.number)))==1) {
		#determination of the condition with the highest number of replicates
		fit.condition<-which(repl.number==max(repl.number))
		if(Verbose) cat("fit.condition",fit.condition,"\n")
	}
	else {
		#more than one condition has the highest number of replicates, therefore the one giving the best fit is chosen
		max.indexes<-which(repl.number==max(repl.number))
		adj.r<-array(,dim=length(max.indexes))
		for(i in 1:length(adj.r)) {
			adj.r[i]<-as.numeric(plgem.fit(data=esdata,fit.condition=max.indexes[i],p=10,q=0.5,verbose=FALSE)$ADJ.R2.MP)
		}
		if(Verbose) cat("adj.r",adj.r,"\n")
		fit.condition<-max.indexes[which(adj.r==max(adj.r))]
		if(length(fit.condition)>1) {
			fit.condition<-fit.condition[1]
			warning("PLGEM fits equally well on more than one condition. \n Condition ",condition.name[fit.condition]," used.\n")
		}
		else {
			if(Verbose) cat("condition ",condition.name[fit.condition]," used \n")
		}
	}

	# fitting and evaluating plgem
	plgemFit<-plgem.fit(data=esdata,covariateNumb,fit.condition,p=10,q=0.5,fittingEval=fitting.eval,plot.file=plotFile,verbose=Verbose)
	# computing observed STN statistics
	obs.stn<-plgem.obsStn(esdata,plgemFit,covariateNumb,verbose=Verbose,baseline.condition=baselineCondition)

	if(repl.number[fit.condition]<3) {
		# since not enough replicates are available for resampling, selection of DEG will be based on ranking
		cat("Less than 3 replicates found in dataset: ranking genes \n")
		DEG.list<-list()
		col.counter<-0
		for(i in (1:condition.number)[-baselineCondition]) {
			col.counter<-col.counter+1
			rankedIDs<-names(sort(abs(obs.stn[,col.counter]),decreasing=TRUE))[1:rank]
			DEG.list[[condition.name[i]]]<-obs.stn[rankedIDs,col.counter]
			names(DEG.list[[condition.name[i]]])<-rankedIDs
		}
	}
	else {
		# computing resampled STN statistics
		res.stn<-plgem.resampledStn(esdata,plgemFit,covariateNumb,iterations=Iterations,baseline.condition=baselineCondition,verbose=Verbose)
    # calculate p-value
    pValues <- plgem.pValue(obs.stn, res.stn, verbose=Verbose)
		# DEG selection
		DEG.list<-plgem.deg(obs.stn,pValues,delta=signLev,verbose=Verbose)
	}

    if(writeFiles) plgem.write.summary(x = DEG.list, verbose = Verbose) # writing DEG list(s) on the disk

	return(DEG.list)
}

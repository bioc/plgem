"run.plgem" <-
function(esdata, signLev=0.001, rank=100, baselineCondition=1, Iterations="automatic", fitting.eval=TRUE, plotFile=FALSE, Verbose=FALSE) {
	library(Biobase)

	#some checks
	if(class(esdata)!="exprSet") stop("Object esdata in function run.plgem is not of class exprSet")
	if(!("conditionName" %in% colnames(pData(esdata)))) stop("the covariate conditionName is not defined in the esdata exprSet")
	if(class(signLev)!="numeric" && class(signLev)!="integer") stop("Argument signLev in function run.plgem is not of class numeric or integer")
	if(class(rank)!="numeric" && class(rank)!="integer") stop("Argument rank in function run.plgem is not of class numeric or integer")
	if(class(baselineCondition)!="numeric" && class(baselineCondition)!="integer") stop("Argument baselineCondition in function run.plgem is not of class numeric or integer")
	if(Iterations!="automatic" && class(Iterations)!="numeric" && class(Iterations)!="integer") stop("Argument Iterations in function run.plgem is neither of class numeric (or integer) nor equal to 'automatic'")
	if(class(fitting.eval)!="logical") stop("Object fitting.eval in function run.plgem is not of class logical")
	if(class(plotFile)!="logical") stop("Object plotFile in function run.plgem is not of class logical")
	if(class(Verbose)!="logical") stop("Argument Verbose in function run.plgem is not of class logical")

	condition.name<-unique(as.character(esdata$conditionName))
	condition.number<-length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object esdata for function run.plgem")
	
	#determining the number of replicates of each condition
	repl.number<-array(,dim=length(condition.name))
	for (i in 1:length(condition.name)) {
	    repl.number[i]<-length(which(esdata$conditionName==condition.name[i]))
	}
	if(max(repl.number)==1) stop ("PLGEM can not be fitted without replicates")
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
	plgemFit<-plgem.fit(data=esdata,fit.condition,p=10,q=0.5,fittingEval=fitting.eval,plot.file=plotFile,verbose=Verbose)
	# computing observed STN statistics
	obs.stn<-plgem.obsStn(esdata,plgemFit,verbose=Verbose,baseline.condition=baselineCondition)

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
		res.stn<-plgem.resampledStn(esdata,plgemFit,iterations=Iterations,baseline.condition=baselineCondition,verbose=Verbose)
		# DEG selection
		DEG.list<-plgem.deg(obs.stn,res.stn,delta=signLev,verbose=Verbose)
	}
	return(DEG.list)
}
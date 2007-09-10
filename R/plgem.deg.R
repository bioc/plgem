"plgem.deg" <-
function(observedStn,plgemResampledStn,delta=0.001,verbose=FALSE) {
	#some checks
	if(class(observedStn)!="matrix") stop("Object observedStn in function plgem.deg is not of class matrix ")
	if(class(plgemResampledStn)!="list") stop("Object plgemResampledStn in function plgem.deg is not of class list ")
	if(class(plgemResampledStn$RESAMPLED.STN)!="matrix") stop("Object plgemResampledStn$RESAMPLED.STN in function plgem.deg is not of class matrix ")
	if(class(delta)!="numeric") stop("Argument delta in function plgem.deg is not of class numeric")
	if(!(all(delta>0) && all(delta<1))) stop("Argument delta in function plgem.deg is not in (0,1)")
	if(class(verbose)!="logical") stop("Argument verbose in function plgem.deg is not of class logical")

	if(verbose) cat("selecting significant DEG:")
	library(Biobase)

	#preparing...
	repl.number<-plgemResampledStn$REPL.NUMBER
	resampledStn<-plgemResampledStn$RESAMPLED.STN
	condition.name<-colnames(observedStn)
	condition.number<-length(condition.name)
	if(verbose) cat("found",condition.number,"condition(s) compared to the baseline.\n")
   	delta.name<-as.character(delta)
	DEG.list<-list()
	up.thr<-array(,dim=c(length(delta),ncol(resampledStn)))
	colnames(up.thr)<-colnames(resampledStn)
	down.thr<-array(,dim=c(length(delta),ncol(resampledStn)))
	colnames(down.thr)<-colnames(resampledStn)

	#calculating up and down thresholds
	for (j in 1:ncol(resampledStn)) {
	    up.thr[,j]<-quantile(resampledStn[,j],1-delta/2,na.rm=TRUE)
	    down.thr[,j]<-quantile(resampledStn[,j],delta/2,na.rm=TRUE)
	}

	#identification of differentially expressed genes (DEG)
	geneIDs <- rownames(observedStn)
	for (i in 1:length(delta)) {
		if(verbose) cat("Delta = ",delta[i],"\n")
		DEG.list[[delta.name[i]]]<-list()
		for (j in 1:condition.number) {
			if(verbose) cat("	Condition = ",condition.name[j],"\n")
			UP<-up.thr[i,as.character(repl.number[condition.name[j]])]
			DOWN<-down.thr[i,as.character(repl.number[condition.name[j]])]
			#selecting DEG
			DEG.index<-which(observedStn[,j]>UP | observedStn[,j]<DOWN)
			DEG.number<-length(DEG.index)
			DEG.stn<-observedStn[DEG.index,j]
			names(DEG.stn)<-geneIDs[DEG.index]
			if(DEG.number==0) {
				DEG.list[[delta.name[i]]][[condition.name[j]]]<-NA
			}
			else {
				DEG.list[[delta.name[i]]][[condition.name[j]]]<-DEG.stn
			}
			if(verbose) {cat("delta:",delta[i],"condition:",condition.name[j],"found",DEG.number,"DEG\n")}
		}
	}
	if(verbose) cat("done with selecting significant DEG.\n\n")
	gc()
	return(DEG.list)
}

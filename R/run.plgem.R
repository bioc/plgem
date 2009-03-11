"run.plgem" <-
function(esdata, signLev=0.001, rank=100, covariate=1, baselineCondition=1, Iterations="automatic", trimAllZeroRows=FALSE, zeroMeanOrSD=c("replace", "trim"), fitting.eval=TRUE, plotFile=FALSE, writeFiles = FALSE, Verbose=FALSE) {

	#some checks
	.checkExpressionSet(esdata)
	if(class(signLev)!="numeric" && class(signLev)!="integer") stop("Argument 'signLev' is not of class 'numeric' or 'integer'")
	if(!(signLev >= 0 && signLev <= 1)) stop("Argument 'signLev' is not in the range [0,1]")
	if(class(rank)!="numeric" && class(rank)!="integer") stop("Argument 'rank' is not of class 'numeric' or 'integer'.")

  covariate <- .checkCovariate(covariate, pData(esdata))
  condition.names <- as.character(pData(esdata)[, covariate])
  baselineCondition <- .checkCondition(baselineCondition, "baselineCondition", covariate, pData(esdata))

	if(Iterations!="automatic" && class(Iterations)!="numeric" && class(Iterations)!="integer") stop("Argument 'Iterations' is neither of class 'numeric' (or 'integer') nor equal to 'automatic'.")
	if(class(fitting.eval)!="logical") stop("Argument 'fitting.eval' is not of class 'logical'.")
	if(class(plotFile)!="logical") stop("Argument 'plotFile' is not of class 'logical'.")
	if(class(Verbose)!="logical") stop("Argument 'Verbose' is not of class 'logical'.")

	condition.name <- unique(condition.names)
	condition.number <- length(condition.name)
	if (condition.number < 2) stop("At least 2 conditions are needed in object 'esdata' for function 'run.plgem'.")

	#determining the number of replicates of each condition
  repl.number <- table(condition.names)
  if(max(repl.number) < 2) stop ("PLGEM can not be fitted without replicates.")

	#determination of the best condition on which to fit the model
	if (Verbose) cat("determining the condition on which to fit the model...\n")
	if (length(which(repl.number == max(repl.number))) == 1) {
		#determination of the condition with the highest number of replicates
		fitCondition <- which(repl.number == max(repl.number))
	} else {
		#more than one condition has the highest number of replicates,
    #therefore the one giving the best fit is chosen
		if(Verbose) cat("more than one condition has the highest number of replicates.\n")
    max.indexes <- which(repl.number == max(repl.number))
    adj.r2 <- sapply(max.indexes, function(x) {
      as.numeric(plgem.fit(data=esdata, covariate=covariate, fitCondition=x,
      p=10, q=0.5, trimAllZeroRows=trimAllZeroRows,
      zeroMeanOrSD=zeroMeanOrSD, verbose=FALSE)$ADJ.R2.MP)
    })
 		if(Verbose) {
      cat("adj. r^2:\n")
      print(adj.r2)
    }
		fitCondition <- max.indexes[which(adj.r2 == max(adj.r2))]
		if(length(fitCondition) > 1) {
			fitCondition <- fitCondition[1]
			warning("PLGEM fits equally well on more than one condition. Condition ",
        sQuote(condition.name[fitCondition]), " used.\n")
		}
	}
	fitCondition <- condition.name[fitCondition]
 	if (Verbose) cat("condition", sQuote(fitCondition), "used.\n")

	# fitting and evaluating plgem
	plgemFit <- plgem.fit(data=esdata, covariate=covariate,
    fitCondition=fitCondition, p=10, q=0.5, trimAllZeroRows=trimAllZeroRows,
    zeroMeanOrSD=zeroMeanOrSD[1], fittingEval=fitting.eval, plot.file=plotFile,
    verbose=Verbose)
	# computing observed STN statistics
	obs.stn<-plgem.obsStn(esdata, plgemFit, covariate, verbose=Verbose, baselineCondition=baselineCondition)

	if(repl.number[fitCondition]<3) {
		# since not enough replicates are available for resampling, selection of DEG will be based on ranking
		cat("Less than 3 replicates found in dataset: ranking genes \n")
		DEG.list<-list()
		col.counter<-0
		for(i in setdiff(condition.name, baselineCondition)) {
			col.counter <- col.counter+1
			DEG.list[[paste("first", rank, sep="")]][[i]] <-
        names(sort(abs(obs.stn[["PLGEM.STN"]][, col.counter]),
        decreasing=TRUE))[1:rank]
		}
		output <- c(obs.stn, list("p.value"=NULL, "significant"=DEG.list))
	}
	else {
		# computing resampled STN statistics
		res.stn <- plgem.resampledStn(esdata, plgemFit, covariate, iterations=Iterations, baselineCondition=baselineCondition, verbose=Verbose)
    # calculate p-value
    pValues <- plgem.pValue(obs.stn, res.stn, verbose=Verbose)
		# DEG selection
		output <- plgem.deg(obs.stn, pValues, delta=signLev, verbose=Verbose)
	}

  if(writeFiles) plgem.write.summary(x = DEG.list, verbose = Verbose) # writing DEG list(s) on the disk

	return(output)
}

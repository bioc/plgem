"plgem.deg" <-
function(observedStn, plgemPval, delta=0.001, verbose=FALSE) {
	#some checks
	if(!is(observedStn, "list")) {
    stop("Object 'observedStn' is not of class 'list'.")
  }
	if(!is(observedStn$PLGEM.STN, "matrix")) {
    stop("Object 'observedStn$PLGEM.STN' is not of class 'matrix'.")
  }
	if(!is(plgemPval, "matrix")) stop("Object 'plgemPval' is not of class 'matrix'.")
	if(!is(delta, "numeric")) stop("Argument 'delta' is not of class 'numeric'.")
	if(!(all(delta>0) && all(delta<1))) stop("One or more elements in argument 'delta' is outside allowed range.")
	if(!is(verbose, "logical")) stop("Argument 'verbose' is not of class 'logical'.")

	if(verbose) cat("selecting significant DEG:")

	#preparing...
	condition.name <- colnames(observedStn$PLGEM.STN)
	condition.number <- length(condition.name)
	if(verbose) cat("found", condition.number, "condition(s) compared to the baseline.\n")
 	delta.name <- as.character(delta)
	DEG.list<-list()

	#identification of differentially expressed genes (DEG)
	geneIDs <- rownames(plgemPval)
	for (i in 1:length(delta)) {
		if(verbose) cat("Delta = ", delta[i], "\n")
		DEG.list[[delta.name[i]]] <- list()
		for (j in 1:condition.number) {
			if(verbose) cat("	Condition = ", condition.name[j], "\n")
			#selecting DEG
			DEG <- geneIDs[plgemPval[, j] < delta[i]]
			
			if(length(DEG) == 0) {
				DEG.list[[delta.name[i]]][[condition.name[j]]]<-NA
			}
			else {
        DEG.list[[delta.name[i]]][[condition.name[j]]]<-DEG
			}
			if(verbose) {cat("delta:", delta[i], "condition:", condition.name[j],
        "found", length(DEG), "DEG\n")}
		}
	}
	if(verbose) cat("done with selecting significant DEG.\n\n")
	gc()
	return(c(observedStn, list("p.value"=plgemPval, "significant"=DEG.list)))
}

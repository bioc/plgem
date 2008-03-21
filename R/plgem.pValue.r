"plgem.pValue" <-
function(observedStn, plgemResampledStn, verbose=FALSE) {

	#some checks
	if(class(observedStn)!="matrix") {
    stop("Object 'observedStn' is not of class 'matrix'.")
  }
	if(class(plgemResampledStn)!="list") {
    stop("Object 'plgemResampledStn' is not of class 'list'.")
  }
	if(class(plgemResampledStn$RESAMPLED.STN)!="matrix") {
    stop("Object 'plgemResampledStn$RESAMPLED.STN' is not of class 'matrix'.")
  }
	if(class(verbose)!="logical") {
    stop("Argument 'verbose' is not of class 'logical'.")
  }

	if(verbose) cat("calculating PLGEM p-values... ")

  # preparing
  repl.number <- plgemResampledStn$REPL.NUMBER
  CDF <- apply(plgemResampledStn$RESAMPLED.STN, 2, ecdf)
  condition.name <- colnames(observedStn)
  pval <- array(, dim=dim(observedStn), dimnames=dimnames(observedStn))
  
  for (i in 1:length(condition.name)) {
    cdf <- CDF[[as.character(repl.number[condition.name[i]])]]
    quant <- cdf(observedStn[, i])
    pval[, i] <- pmin(2 * quant, 2 * (1 - quant))
  }
  
  if(verbose) cat("done.\n\n")
  return(pval)

}

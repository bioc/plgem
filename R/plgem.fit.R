"plgem.fit" <- function(data, covariate=1, fitCondition=1, p=10, q=0.5,
  trimAllZeroRows=FALSE, zeroMeanOrSD=c("replace", "trim"), fittingEval=FALSE,
  plot.file=FALSE, prefix=NULL, gPar=setGpar(), verbose=FALSE) {
  
	# some checks..
	.checkExpressionSet(data)
  covariate <- .checkCovariate(covariate, pData(data))
  fitCondition <- .checkCondition(fitCondition, "fitCondition", covariate,
    pData(data))
  fittingEvalFileName <- paste(as.character(prefix), "fittingEval.png", sep="_")
  .checkExistence(fittingEvalFileName)

	if(!is(p, "numeric") && !is(p, "integer")) {
    stop("Argument 'p' is not of class 'numeric' or 'integer'")
  }
	if(!is(q, "numeric") && !is(q, "integer")) {
    stop("Argument 'q' is not of class 'numeric' or 'integer'")
  }
	if(!(q >= 0 && q <= 1)) stop("Argument 'q' is not in the range [0, 1]")
	if(!is(fittingEval, "logical")) {
    stop("Argument 'fittingEval' is not of class 'logical'")
  }
	if(!is(plot.file, "logical")) {
    stop("Argument 'plot.file' is not of class 'logical'")
  }
	if(!is(verbose, "logical")) {
    stop("Argument 'verbose' is not of class 'logical'")
  }
	if(verbose) cat("Fitting PLGEM...\n")

	condition.names <- as.character(pData(data)[, covariate])
	if(verbose) cat("samples extracted for fitting:\n")
	data <- data[, condition.names == fitCondition]
	if(verbose) print(pData(data))
	if(length(sampleNames(data)) < 2) {
    stop("At least 2 replicates needed to fit PLGEM")
  }
  
	# extract data matrix and optionally trim rows containing only zero values
	dataMatrix <- exprs(data)
  if (trimAllZeroRows) {
    allZero <- apply(dataMatrix, 1, function(z) all(z==0) )
    if (sum(allZero)) {
      if (verbose) cat("trimming", sum(allZero),
        "rows with only zero values.\n")
      dataMatrix <- dataMatrix[!allZero, ]
    }
  }
	
	# calculate row means and standard deviations
	stats <- as.data.frame(t(apply(dataMatrix, 1, function(x) {
    c("mean"=mean(x, na.rm=T), "sd"=sd(x, na.rm=T))
  })))
  stats <- stats[order(stats$mean), ]
  
  # decide what to do with rows with non-positive mean or zero standard deviation
  zeroMeanOrSD <- match.arg(zeroMeanOrSD)
  switch(zeroMeanOrSD,
    trim = {
      toRemove <- union(which(stats$mean <= 0), which(stats$sd == 0))
      if (length(toRemove) > 0) {
        if (verbose) cat("removing", length(toRemove),
          "rows with non-positive mean or zero standard deviation...\n")
        stats <- stats[-toRemove, ]  
      }
    },
    replace = {
      toRemove <- which(stats$mean <= 0)
      if (length(toRemove) > 0) {
        if (verbose) cat("replacing", length(toRemove),
          "non-positive means with smallest positive mean...\n")
        stats$mean <- replace(stats$mean, stats$mean <= 0,
          min(stats$mean[stats$mean > 0]))
      }
      toRemove <- which(stats$sd == 0)
      if (length(toRemove) > 0) {
        if (verbose) cat("replacing", length(toRemove),
          "zero standard deviations with smallest non-zero standard deviation...\n")
        stats$sd <- replace(stats$sd, stats$sd == 0, min(stats$sd[stats$sd > 0]))
      }
    }  
  )

	# determine modelling points (MP)
	if(verbose) cat("determining modelling points...\n")
	MP.x <- array(, dim=p)
	MP.y <- array(, dim=p)
	limits <- round(seq(1, nrow(stats), length.out=p + 1))
	for (i in 1:p) {
		MP.y[i] <- quantile(stats$sd[limits[i]:limits[i + 1]], q, na.rm=TRUE)	
		MP.x[i] <- median(stats$mean[limits[i]:limits[i + 1]], na.rm=TRUE)
	}

	# fit linear regression over the modelling points
	if(verbose) cat("fitting data and modelling points...\n")
	MP.lm <- lm(log(as.numeric(MP.y)) ~ log(as.numeric(MP.x)))
	slope <- as.numeric(coef(MP.lm)[2])
	if (slope > 1) warning("PLGEM slope is higher than 1")
	if (slope < 0.5) warning("PLGEM slope is lower than 0.5")
	intercept <- as.numeric(coef(MP.lm)[1])
	adj.r2.mp <- summary.lm(MP.lm)$adj.r.squared
	if (adj.r2.mp < 0.95) warning("Adjusted r^2 is lower than 0.95")
	data.pearson <- cor(log(stats$mean), log(stats$sd))
	if (data.pearson < 0.85)
    warning("Pearson correlation coefficient is lower than 0.85")

	if(fittingEval) {
		if(plot.file) png(filename=fittingEvalFileName, width=600, height=600)

		# computation of model residuals
		modeled.spread <- (stats$mean^slope) * (exp(intercept))
		residuals <- (log(stats$sd) - log(modeled.spread))
		maxResidual <- ceiling(max(abs(residuals)))

		layout(matrix(1:4, 2, 2))
		par(mar=c(3, 3, 3.5, 1), mgp=c(1.5, .5, 0))

		# contour plot of data & model fit
		meanSd.kde2d <- kde2d(log(stats$mean), log(stats$sd), n=50)
		if (is.null(gPar$minLnM)) gPar$minLnM <- min(meanSd.kde2d$x)
 		if (is.null(gPar$maxLnM)) gPar$maxLnM <- max(meanSd.kde2d$x)
		if (is.null(gPar$minLnSD)) gPar$minLnSD <- min(meanSd.kde2d$y)
 		if (is.null(gPar$maxLnSD)) gPar$maxLnSD <- max(meanSd.kde2d$y)
		contour(meanSd.kde2d, col=colors()[40:55], nlevels=15, xlab="ln(mean)",
      ylab="ln(sd)", cex=0.5, main="Power Law Global Error Model",
      xlim=c(gPar$minLnM, gPar$maxLnM), ylim=c(gPar$minLnSD, gPar$maxLnSD))
		points(log(MP.x), log(MP.y), cex=2, col="black")
		abline(intercept, slope, col="red", lwd=2)
		
    xy <- par("usr")
		text(xy[1]+0.2*diff(xy[1:2]), xy[3]+0.9*diff(xy[3:4]), paste("slope =",
      signif(slope, 3), "\nintercept =", signif(intercept, 3)), cex=0.8)
		text(xy[1]+0.8*diff(xy[1:2]), xy[3]+0.1*diff(xy[3:4]), paste("adj. r2 =",
      signif(adj.r2.mp, 3), "\nPearson =", signif(data.pearson, 3)), cex=0.8)

		# contour plot of residuals vs rank of mean
		residuals.kde2d <- kde2d(rank(stats$mean, ties.method="first"), residuals, n=50)
		if (is.null(gPar$maxRes)) gPar$maxRes <- maxResidual
 		if (is.null(gPar$minRes)) gPar$minRes <- -maxResidual
		contour(residuals.kde2d, col=colors()[40:55], nlevels=15,
      xlab="Rank of mean", ylab="Residuals", ylim=c(gPar$minRes, gPar$maxRes),
      main="Residuals vs. rank", cex=0.5)

		# histogram of residuals
		hist(residuals, breaks=100, xlab="Residuals", ylab="Counts",
      xlim=c(gPar$minRes, gPar$maxRes))
		box()

		# qqplot of residuals distribution versus standard normal distribution
		qqnorm(residuals, xlab="Standard Normal", ylab="Residuals", cex=0.5)
		
		# add a title to the outer margin of the plot
		mtext(paste("PLGEM fitted on condition '", fitCondition, "'", sep=""),
      line=-1.1, outer=TRUE)
		if(plot.file) dev.off()
	}

	# return model parameters
	if(verbose) cat("done with fitting PLGEM.\n\n")
	gc()
	return(list(SLOPE=slope, INTERCEPT=intercept, DATA.PEARSON=data.pearson,
    ADJ.R2.MP=adj.r2.mp, COVARIATE=covariate, FIT.CONDITION=fitCondition))
}

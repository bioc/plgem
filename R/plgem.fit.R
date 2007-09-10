"plgem.fit" <-
function(data, covariateNumb=1, fit.condition=1,p=10,q=0.5,fittingEval=FALSE,plot.file=FALSE,verbose=FALSE) {
	library(Biobase)
	library(MASS)

	# some checks..
	if(class(data)!="ExpressionSet") stop("Object data in function plgem.fit is not of class ExpressionSet")
	if(covariateNumb > ncol(pData(data))) stop("covariateNumb is greater than the number of covariates in phenodata of data")
	if(class(fit.condition)!="numeric" && class(fit.condition)!="integer") stop("Argument fit.condition in function plgem.fit is not of class numeric or integer")
	if(class(p)!="numeric" && class(p)!="integer") stop("Argument p in function plgem.fit is not of class numeric or integer")
	if(class(q)!="numeric" && class(q)!="integer") stop("Argument q in function plgem.fit is not of class numeric or integer")
	if(class(fittingEval)!="logical") stop("Argument fittingEval in function plgem.fit is not of class logical")
	if(class(plot.file)!="logical") stop("Argument plot.file in function plgem.fit is not of class logical")
	if(class(verbose)!="logical") stop("Argument verbose in function plgem.fit is not of class logical")

	if(verbose) cat("fitting PLGEM..","\n")

	if(verbose) cat("samples extracted for fitting:","\n")
	condition.names<-as.character(pData(data)[,covariateNumb])
	condition.name<-unique(condition.names)
	data<-data[,which(condition.names==condition.name[fit.condition])]
	if(verbose) print(pData(data))
	row<-length(featureNames(data))
	if(length(sampleNames(data))<2) stop("At least 2 replicates needed to fit PLGEM")

	# 'data' mean and standard deviation
	dataMatrix<-exprs(data)
	data.mean<-rowMeans(dataMatrix,na.rm=TRUE)
	data.mean<-sort(data.mean)
	data.sd<-sd(t(dataMatrix),na.rm=TRUE)
	data.sd<-data.sd[names(data.mean)]
	data.mean<-replace(data.mean,data.mean<=0,min(data.mean[data.mean>0]))
	data.sd<-replace(data.sd,data.sd<=0,min(data.sd[data.sd>0]))

	# Determination of Modelling Points (MP)
	if(verbose) cat("determining modelling points...\n")
	limits<-array(,dim=p+1)
	MP.x<-array(,dim=p)
	MP.y<-array(,dim=p)
	for (i in 0:p) {limits[i+1]<-round(row*i/p)}
	for (i in 1:(length(limits)-1)) {
		MP.y[i]<-quantile(data.sd[limits[i]:limits[i+1]],q,na.rm=TRUE)	
		MP.x[i]<-median(data.mean[limits[i]:limits[i+1]],na.rm=TRUE)
	}

	# fit of linear regression over the Modelling Points
	if(verbose) cat("fitting data and modelling points...\n")
	MP.lm<-lm(log(as.numeric(MP.y))~log(as.numeric(MP.x)))
	slope<-coef(MP.lm)[2]
	intercept<-coef(MP.lm)[1]
	adj.r2.mp<-summary.lm(MP.lm)$adj.r.squared
	data.pearson<-cor(log(data.mean),log(data.sd))

	if(fittingEval) {
		if(plot.file) {png(file="fittingEval.png",width=600,height=600)}
		# computation of model-residuals
		modeled.spread<-(data.mean^slope)*(exp(intercept))
		residuals<-(log(data.sd)-log(modeled.spread))
		maxResidual<-ceiling(max(abs(residuals)))

		layout(matrix(1:4, 2, 2))
		# contour plot of data & model fit
		mainTitle<-paste("PLGEM parameters:\nslope = ",signif(slope,3),"\nintercept = ",signif(intercept,3),sep="")
		meanSd.kde2d<-kde2d(log(data.mean), log(data.sd),n=50)
		contour(meanSd.kde2d,col=colors()[40:55],nlevels=15,xlab="ln(mean)",ylab="ln(sd)",cex=0.5,main=mainTitle,cex.main=0.8)
		points(log(MP.x),log(MP.y),cex=2,col="black")
		abline(intercept,slope,col="red",lwd=2)
		# contour plot of residuals vs rank of mean
		residuals.kde2d<-kde2d(rank(data.mean,ties.method="first"),residuals,n=50)
		contour(residuals.kde2d,col=colors()[40:55],nlevels=15,xlab="Rank of mean",ylab="Residuals",ylim=c(-maxResidual,maxResidual),cex=0.5)
		# histogram of residuals
		hist(residuals,breaks=100,xlab="Residuals",ylab="Counts",xlim=c(-maxResidual,maxResidual))
		# qqplot of residuals distribution versus standard normal distribution
		qqnorm(residuals,xlab="Standard Normal",ylab="Residuals",cex=0.5)
		if(plot.file) {dev.off()}
	}

	# return model parameters
	if(verbose) cat("done with fitting PLGEM.\n\n")
	gc()
	return(list(SLOPE=slope,INTERCEPT=intercept,DATA.PEARSON=data.pearson,ADJ.R2.MP=adj.r2.mp,FIT.CONDITION=fit.condition))
}

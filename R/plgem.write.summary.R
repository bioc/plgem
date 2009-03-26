"plgem.write.summary" <- function(x, prefix=NULL, verbose=FALSE) {
  if(class(x)!="list" || !identical(names(x), c("fit", "PLGEM.STN", "p.value", "significant"))) stop("x has be the output of either the plgem.deg or the run.plgem function")
  if (verbose) cat("Writing files...\n")
  fitFileName <- paste(prefix, "fit.csv", sep="_")
  .checkExistence(fitFileName)
  if (verbose) cat("   ", fitFileName, "\n")
  write.table(t(as.data.frame(x$fit)), file=fitFileName, sep=",",
    col.names=FALSE)
  
  if (is.null(x$p.value)) {
    mat <- cbind(x$PLGEM.STN)
    colnames(mat) <- paste(colnames(x$PLGEM.STN), "STN", sep="_")
    stnFileName <- paste(prefix, "stn.csv", sep="_")
    .checkExistence(stnFileName)
    if (verbose) cat("   ", stnFileName, "\n")
    write.csv(mat, file=stnFileName)
  } else {
    mat <- cbind(x$PLGEM.STN, x$p.value)
    colnames(mat) <- c(paste(colnames(x$PLGEM.STN), "STN", sep="_"),
      paste(colnames(x$p.value), "p-value", sep="_"))
    stnFileName <- paste(prefix, "stn_p-value.csv", sep="_")
    .checkExistence(stnFileName)
    if (verbose) cat("   ", stnFileName, "\n")
    write.csv(mat, file=stnFileName)
  }
  
  for (i in names(x$significant)) {
    for (j in names(x$significant[[i]])) {
      fname <- paste(prefix, "_", j, "_", i, ".txt", sep="")
      .checkExistence(fname)
      if (verbose) cat("   ", fname, "\n")
      write(x$significant[[i]][[j]], file=fname)
    }
  }
  
  if (verbose) cat("...to folder", getwd(), "\n")
}

"plgem.write.summary" <-
function(x, prefix=NULL, verbose=FALSE) {
  if(class(x)!="list" || !identical(names(x), c("fit", "PLGEM.STN", "p.value", "significant"))) stop("x has be the output of either the plgem.deg or the run.plgem function")
  if (verbose) cat("Writing files\n")
  write.table(t(as.data.frame(x$fit)), file=paste(prefix, "fit.csv", sep="_"),
    sep=",", col.names=FALSE)
  
  if (is.null(x$p.value)) {
    mat <- cbind(x$PLGEM.STN)
    colnames(mat) <- paste(colnames(x$PLGEM.STN), "STN", sep="_")
    write.csv(mat, file=paste(prefix, "stn.csv", sep="_"))
  } else {
    mat <- cbind(x$PLGEM.STN, x$p.value)
    colnames(mat) <- c(paste(colnames(x$PLGEM.STN), "STN", sep="_"), paste(colnames(x$p.value), "p-value", sep="_"))
    write.csv(mat, file=paste(prefix, "stn_p-value.csv", sep="_"))
  }
  
  f1 <- names(x$significant)
  for (i in f1) {
    y <- x$significant[[i]]
    f2 <- names(y)
    for (j in f2) {
      z <- x$significant[[i]][[j]]
      fname <- paste(prefix, "_", j, "-", i, ".txt", sep="")
      if (verbose) cat(fname, "\n")
      write(z, file=fname)
      }
    }
  if (verbose) cat("to folder", getwd(), "\n")
}

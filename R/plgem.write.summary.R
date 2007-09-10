"plgem.write.summary" <-
function(x, verbose=FALSE) {
  if(class(x)!="list" || class(x[[1]])!="list") stop("x has be the output of either the plgem.deg or the run.plgem functions, i.e. a list of list(s) of named vectors")
  if (verbose) cat("Writing files\n")
  f1 <- names(x)
  for (i in f1) {
    y <- x[[i]]
    f2 <- names(y)
    for (j in f2) {
      z <- x[[i]][[j]]
      fname <- paste(j, "-", i, ".txt", sep="")
      if (verbose) cat(fname, "\n")
      write.table(round(z,3), file=fname, sep="\t", quote=FALSE, col.names=FALSE)
      }
    }
  if (verbose) cat("to folder", getwd(), "\n")
}

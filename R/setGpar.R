setGpar <- function(minLnM=NULL, maxLnM=NULL, minLnSD=NULL, maxLnSD=NULL,
  minRes=NULL, maxRes=NULL) {
  
  if(!missing(minLnM) && class(minLnM) != "numeric" && class(minLnM) != "integer") {
    stop("Argument 'minLnM' is not of class 'numeric' or 'integer'")
  }
  if(!missing(maxLnM) && class(maxLnM) != "numeric" && class(maxLnM) != "integer") {
    stop("Argument 'maxLnM' is not of class 'numeric' or 'integer'")
  }


  list(minLnM=minLnM, maxLnM=maxLnM, minLnSD=minLnSD, maxLnSD=maxLnSD,
  minRes=minRes, maxRes=maxRes)
}

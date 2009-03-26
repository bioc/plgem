# ==========================================================================
# private functions of the 'plgem' package
# ==========================================================================

# check input ExpressionSet
.checkExpressionSet <- function(eset) {
  funCall <- as.character(sys.call(-1))[1]
  if (class(eset) != "ExpressionSet") {
    stop("Input dataset for function ", sQuote(funCall), " is not of class ",
      sQuote("ExpressionSet"), ".")
  }
  if (ncol(pData(eset)) < 1) {
    stop("No covariates defined in the input ExpressionSet for function ",
      sQuote(funCall), ".")
  }
}

# check input covariate for consistency with phenoData
.checkCovariate <- function(covar, pdat) {
  funCall <- as.character(sys.call(-1))[1]
  if (length(covar) > 1) {
    covar <- covar[1]
    warning("Multiple covariates specified in ", sQuote(funCall),
      ". Only the first one will be used.")
  }
  if (class(covar) == "numeric" || class(covar) == "integer") {
    if(covar < 1) stop("Argument ", sQuote("covariate"), " for function ",
      sQuote(funCall), " must be >= 1.")
    if(as.integer(covar) > ncol(pdat)) stop("Argument ", sQuote("covariate"),
      " is greater than the number of covariates in the input ExpressionSet for function ",
      sQuote(funCall), ".")
    return(colnames(pdat)[as.integer(covar)])
  }
  if (class(covar) == "character") {
    if(!(covar %in% colnames(pdat))) stop ("covariate ", sQuote(covar),
      " is not defined in the input ExpressionSet for function ",
      sQuote(funCall), ".")
    return(covar)
  }
  if (class(covar) != "numeric" && class(covar) != "integer" &&
    class(covar) != "character") {
    stop("Argument ", sQuote("covariate"), " for function ", sQuote(funCall),
      " must be one of class ", sQuote("numeric"), ", ", sQuote("integer"),
      " or ", sQuote("character"), ".")
  }
}

# check input condition for consistency with phenoData
.checkCondition <- function(cond, argName, covar, pdat) {
  funCall <- as.character(sys.call(-1))[1]
  vars <- unique(as.character(pdat[, covar]))
  if (length(cond) > 1) {
    cond <- cond[1]
    warning("Multiple conditions specified in ", sQuote(funCall),
      ". Only the first one will be used.")
  }
  if (class(cond) == "numeric" || class(cond) == "integer") {
    if(cond < 1) stop("Argument ", sQuote(argName), " for function ",
      sQuote(funCall), " must be >= 1.")
    if(as.integer(cond) > length(unique(vars))) {
      stop("Argument ", sQuote(argName),
        " is greater than the number of conditions in the input ExpressionSet for function ",
      sQuote(funCall), ".")
    }
    return(vars[as.integer(cond)])
  }
  if (class(cond) == "character") {
    if(!(cond %in% vars)) {
      stop("condition ", sQuote(cond),
      " is not defined in the input ExpressionSet for function ",
        sQuote(funCall), ".")
    }
    return(cond)
  }
  if (class(cond) != "numeric" && class(cond) != "integer" && class(cond) != "character") {
    stop("Argument ", sQuote(argName), " for function ", sQuote(funCall),
      " must be one of class ", sQuote("numeric"), ", ", sQuote("integer"),
      " or ", sQuote("character"), ".")
  }
}

# get signal-to-noise ratio from location and spread
.stn <- function(location1, location2, spread1, spread2) {
	(location2 - location1)/(spread1 + spread2)
}

# get PLGEM-based estimate of spread, based on location
.plgemSpread <- function(location, slope, intercept) {
	(location^slope) * exp(intercept)
}

# check for file existence before saving
.checkExistence <- function(fileName) {
  funCall <- as.character(sys.call(-1))[1]
  if (file.exists(fileName)) {
    choice <- menu(c("Overwrite?", "Abort?"), graphics=TRUE,
      title=paste("File", fileName, "already exists"))
    if (choice == 2) stop(sQuote(funCall), " aborted.")
  }
}


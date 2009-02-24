# ==========================================================================
# private functions of the 'plgem' package
# ==========================================================================

# check input covariate for consistency with phenoData
.checkCovariate <- function(covar, pdat) {
  if (length(covar) > 1) {
    covar <- covar[1]
    warning("Multiple covariates specified. Only the first one will be used.")
  }
  if (class(covar) == "numeric" || class(covar) == "integer") {
    if(covar < 1) stop("Argument 'covariate' must be >= 1.")
    if(as.integer(covar) > ncol(pdat)) stop("Argument 'covariate' is greater than the number of covariates in the ExpressionSet.")
    return(colnames(pdat)[as.integer(covar)])
  }
  if (class(covar) == "character") {
    if(!(covar %in% colnames(pdat))) stop ("covariate '", covar, "' is not defined in the ExpressionSet.")
    return(covar)
  }
  if (class(covar) != "numeric" && class(covar) != "integer" && class(covar) != "character") {
    stop("Argument 'covariate' must be one of class 'numeric', 'integer' or 'character'.")
  }
}

# check input condition for consistency with phenoData
.checkCondition <- function(cond, argName, vars) {
  if (length(cond) > 1) {
    cond <- cond[1]
    warning("Multiple conditions specified. Only the first one will be used.")
  }
  if (class(cond) == "numeric" || class(cond) == "integer") {
    if(cond < 1) stop("Argument '", argName, "' must be >= 1.")
    if(as.integer(cond) > length(unique(vars))) stop("Argument '", argName, "' is greater than the number of covariates in the ExpressionSet.")
    return(unique(vars)[as.integer(cond)])
  }
  if (class(cond) == "character") {
    if(!(cond %in% vars)) stop("condition '", cond, "' is not defined in the ExpressionSet.")
    return(cond)
  }
  if (class(cond) != "numeric" && class(cond) != "integer" && class(cond) != "character") {
    stop("Argument '", argName, "' must be one of class 'numeric', 'integer' or 'character'.")
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


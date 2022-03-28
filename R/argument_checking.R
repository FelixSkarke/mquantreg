arg_check <- function(q,
                      k,
                      method,
                      scale.estimator,
                      maxit,
                      X,
                      offset,
                      case.weights,
                      var.weights,
                      weights.x,
                      qgrid) {


##############################################
## Checking arguments ##
## Check if valid m-quantile q is chosen
    if (any(q <= 0) | any (q >= 1)){
      stop("Choosen M-quantile is not valid. q must be greater than 0 and less than 1.")
      }


##############################################
## Check if valid method is selected
   if (!(method %in% c("poisson", "quasipoisson", "nb", "binom", "continuous")))
      {stop("Choosen method is not implemented. Please select a valid method.
          Avaiable methods are: poisson, quasipoisson, nb, binom, continuous"
           )
      }

##############################################
## Check if design matrix is singular
    if (qr(X)$rank < ncol(X)){
      stop("X matrix is singular")
    }

##############################################
## Check offset is valid


# Check for correct offset length
    if (length(offset) != nrow(X)){
      stop("Length of offset must equal the number of observations")
    }


##############################################
## Check if passed case-weights are valid

    if (length(case.weights) != nrow(X)){
      stop("Length of case weights must equal the number of observations")
      }

    if (any(case.weights < 0)){
      stop("Negative weights are not allowed")
      }

##############################################
## Check if maxit is valid
    if ((maxit < 0) | (maxit %% 1 != 0)){
      stop("Number of iterations must be a positive integer.")
    }


##############################################
## Check if k is valid

    if (k <= 0) {
      stop("k must be a number larger than 0.")
    }

    if ((k == Inf) && (method %in% c("quasipoisson", "nb"))) {
      stop("Expectile estimation not implemented. Choose a finite k.")
      }

##############################################
## Check if passed case-weights are valid

    if (length(var.weights) != nrow(X)) {
      stop("Length of var.weights must equal number of observations.")
      }

    if (any(var.weights < 0)) {
      stop("Negative var.weights are not allowed")
      }


##############################################
## Check if implemented scale estimator passed
    if (!((scale.estimator == "Mad") || (scale.estimator == "cMad"))) {
      stop("The chosen scale estimator is not implemented. Please choose Mad or cMad.")
    }

##############################################
## Check if passed weights.x are valid
    if (!isTRUE(weights.x) && !isFALSE(weights.x)) {
      stop("The argument weight.x has to be either TRUE or FALSE.")
    }


#############################################
## Check if passed qvalues for qscores are valid
    if (any(qgrid < 0) || any(qgrid > 1)) {
      stop("Q-values for computing qscores must be between 0 and 1.")
    }
  }





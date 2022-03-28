#' M-Quantile Regression
#'
#'The mquantreg-function implements m-quantile regression for multiple types of
#'dependent data.It serves as an interface encompassing all types of m-quantile
#'regression implemented in this package. It is recommended to always use the mquantreg-
#'function, although it is possible to call the models from their respective functions.
#'
#'@param formula a formula object, with the response on the left of a \code{~} operator
#'and the covariates, separated by \code{+} operators on the right hand side.
#'@param data a data.frame, which contains the variables named in the formula
#'and in the subset and the weights arguments.
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'Can also be a vector of m-quantiles.
#'@param k a number greater than 0. k is the parameter for the huber psi-function and
#'defaults to 1.345.
#'@param method string, which specifies the model to estimate.
#'The available models are: \cr
#'\itemize{
#'\item \code{"continuous"} - Continuous M-Quantile Regression. See \link{mqcont}.
#'\item \code{"binom"} - Binary M-Quantile Regession. See \link{mqbinom}.
#'\item \code{"poisson"} - Poisson M-Quantile Regession for Count Data. See \link{mqpoisson}.
#'\item \code{"quasipoisson"} - Quasipoisson M-Quantile Regession for Count Data. See \link{mqquasipoisson}.
#'\item \code{"nb"} - Negative Binomial  M-Quantile Regression for Count Data. See \link{mqnb}.
#'}
#'@param scale.estimator string, defines which function is used to estimate the
#'scale parameter used for peliminary estimation in the continuous regression.
#'The following scale estimators are available: \cr
#'\itemize{
#'\item \code{"Mad"} - Mean Absolute Deviation (default)
#'\item \code{"cMad"} - Corrected Mean Absolute Deviation
#'}
#'@param na.action defaults to \code{na.fail}.
#'@param maxit a integer, which defines the maximum number of iterations
#'before the estimation process is halted. If the maximum number of iterations is
#'reached a warning message will be prompted.
#'@param offset a numeric vector of length equal to the number of observations.
#'This can be used to specify an a priori known component to be
#'included in the linear predictor during fitting. Defaults to \code{NULL}.
#'@param case.weights a numeric vector of observation weights;
#'if supplied, the algorithm fits to minimize the sum of the weights multiplied into
#'the absolute residuals.The length of case.weights must be the same as the number of
#'observations. The weights must be nonnegative and it is recommended that they be
#'strictly positive.Defaults to \code{NULL}.
#'@param weights.x boolean, for the estimation of robust glm to
#'obtain initial coefficiants in the discrete models.Defaults to \code{FALSE}.
#'@param var.weights dataframe, if supplied the residuals are scaled. Has to be of
#'same length as number of observations.
#'@param acc defines convergence criteria.
#'@param compute.qscores boolean, controls if qscores are to be computed.Defaults to
#'\code{FALSE}.
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@return mquantreg returns an object of class "mquantfit". The class has a number of avaiable methods:\cr
#'\code{summary()}, \code{print()}, \code{fitted()}, \code{qscore()}, \code{huberise()} and \code{r2pseudo() }.
#'An object of class "mquantfit" is a list containing at least the following components:
#'\itemize{
#'\item coefficients - a named vector of coefficients.
#'\item residuals - the residuals, that is response minus fitted values.
#'\item fitted.values - the fitted mean values.
#'\item q.scores - a vector with the q-score of each observation
#'\item var.beta - variance of the estimators.
#'\item scale - estimated scale parameter.
#'}
#'Additional components are available, including all arguments passed when calling the mquantreg-function.
#' @references
#'\itemize{
#'\item Breckling, J., & Chambers, R. (1988). M-quantiles. Biometrika, 75(4), 761-771.
#'\item Bianchi, A., Fabrizi, E., Salvati, N., & Tzavidis, N. (2015). M-quantile Regression:
#'Diagnostics and the parametric representation of the model.Book of Abstracts SIS - CLADAG, 303-306
#'}
#'
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'
#'                   measurement.error = 0.5)
#'fit <- mquantreg(formula = Y ~ x1 + x2, data = df, q  = 0.05, method = "continuous")
#'summary(fit)
#'@export
#'@import MASS ggplot2
mquantreg <- function(formula,
                      data,
                      q = 0.5,
                      k = 1.345,
                      method = "continuous",
                      scale.estimator = "Mad",
                      na.action = na.fail,
                      maxit = 1000,
                      offset = NULL,
                      case.weights = NULL,
                      var.weights = NULL,
                      weights.x = FALSE,
                      acc = 1e-04,
                      compute.qscores = FALSE,
                      qgrid = c(0.1,
                                0.15,
                                0.2,
                                0.25,
                                0.30,
                                0.35,
                                0.4,
                                seq(from = 0.45, to = 0.55, by = 0.005),
                                0.60,
                                0.65,
                                0.70,
                                0.75,
                                0.8,
                                0.85,
                                0.9)
                      ) {

  ##############################################
  ## Creating X and Y from formula ##
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(x = c("formula", "data", "weights", "na.action"),
              table = names(mf),
              nomatch = 0L
  )

  mf                    <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]]              <- quote(stats::model.frame)
  mf                    <- eval(expr = mf, envir = parent.frame())


  mt <- attr(mf, "terms")
  Y  <- model.response(mf, "any")

  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y)<- NULL

    if (!is.null(nm)){
      names(Y) <- nm
    }
  }


  X <- if (!is.empty.model(mt)){
    model.matrix(mt, mf)
  }else{
    matrix(data = NA, nrow = NROW(Y), 0)
  }

  # preparing weights and offset:
  # Check for models with additiv offset
  if (is.null(offset) && (method %in% c("binom", "continuous"))){
    offset <- rep(0, nrow(X))
  }

  # Check for models with multiplicative offset
  if (is.null(offset) && (method %in% c("poisson", "quasipoisson", "nb"))){
    offset <- rep(1, nrow(X))
  }


  if (is.null(case.weights)){
    case.weights <- rep(1, nrow(X))
  }


  if (is.null(var.weights)) {
    var.weights <- rep(1, nrow(X))
  }

    # argument checking
    arg_check(q               = q,
              method          = method,
              X               = X,
              offset          = offset,
              case.weights    = case.weights,
              maxit           = maxit,
              k               = k,
              var.weights     = var.weights,
              scale.estimator = scale.estimator,
              weights.x       = weights.x,
              qgrid           = qgrid)

  ##############################################
  ## Run estimation according to chosen method
  fit <- switch(method,
                poisson = mqpoisson(x = X,
                                    y = Y,
                                    offset = offset,
                                    case.weights = case.weights,
                                    maxit = maxit,
                                    acc = 1e-04,
                                    var.weights = var.weights,
                                    weights.x = weights.x,
                                    q = q,
                                    k = k,
                                    compute.qscores = compute.qscores
                                    ),
                quasipoisson = mqquasipoisson(x = X,
                                              y = Y,
                                              offset = offset,
                                              case.weights = case.weights,
                                              maxit = maxit,
                                              acc = 1e-04,
                                              var.weights = var.weights,
                                              weights.x = weights.x,
                                              q = q,
                                              k = k,
                                              phi.init = 1,
                                              compute.qscores = compute.qscores
                                              ),
                  nb = mqnb(x = X,
                            y = Y,
                            offset = offset,
                            case.weights = case.weights,
                            maxit = maxit,
                            acc = 1e-04,
                            var.weights = var.weights,
                            weights.x = weights.x,
                            q = q,
                            k = k,
                            theta.init = 1,
                            compute.qscores = compute.qscores
                            ),
                  binom = mqbinom(x = X,
                                  y = Y,
                                  case.weights = case.weights,
                                  maxit = maxit,
                                  acc = 1e-04,
                                  var.weights = var.weights,
                                  weights.x = weights.x,
                                  q = q,
                                  k = k,
                                  compute.qscores = compute.qscores
                                  ),
                  continuous = mqcont(x = X,
                                      y = Y,
                                      offset = offset,
                                      case.weights = case.weights,
                                      maxit = maxit,
                                      acc = 1e-04,
                                      var.weights = var.weights,
                                      q = q,
                                      k = k,
                                      compute.qscores = compute.qscores,
                                      scale.estimator = scale.estimator,
                                      qgrid = qgrid)

  )

  # ##############################################
  ## Define output
  # Create mquantfit class
  class(fit) <- "mquantfit"

  fit$data                   <- data
  fit$X                      <- X
  fit$Y                      <- Y
  fit$method                 <- method
  fit$formula                <- formula
  fit$terms                  <- mt
  fit$xlevels                <- .getXlevels(mt, mf)
  fit$call                   <- cl
  fit$case.weights           <- case.weights
  fit$var.weights            <- var.weights
  fit$weights.x              <- weights.x
  fit$residuals              <- fit$residuals
  fit$k                      <- k
  fit$fitted.values          <- fit$fitted.values
  fit$scale.estimator        <- scale.estimator
  fit$coefficients           <- cbind(q, t(fit$coefficients))
  colnames(fit$coefficients) <- c("q", colnames(X))
  fit$compute.qscores        <- compute.qscores

  return(fit)
}


#' Print M-Quantile Estimates
#'
#' The Print-method returns the coefficients for all estimated m-quantile in one table.
#'
#'@param x a mquantfit object, produced by calling the mquantreg-function.
#'@param ... further arguments passed to or from other methods.
#'\code{summary()}, \code{print()}, \code{fit()} and \code{predict()}-methods are avaiable.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'print(fit)
#'@export
print.mquantfit <- function(x, ...){
    cat("Call:\n")
    dput(x$call)

    coef = coef(x)
    cat("\nCoefficients:\n")
    print(coef)
    rank = x$rank
    nobs = nrow(residuals(x))
    if (is.matrix(coef))
      p = dim(coef)[2] - 1
    else
      p = length(coef)
    rdf = nobs - p
    cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
    if (!is.null(attr(x, "na.message")))
      cat(attr(x, "na.message"), "\n")
    invisible(x)
  }


#' Regression Output for Estimated M-Quantile Model
#'
#' The Summary()-method prints the regression output. The regression output includes
#' the coefficients and standard errors. Furthermore, information about the residuals and
#' the Huberisation are printed. For continuous estimation also pseudo R-square is computed.
#' The regression output for each m-quantile is returned seperately.
#'
#'@param object a mquantfit object, produced by calling the mquantreg-function.
#'@param ... additional arguments affecting the summary produced
#'\code{summary()}, \code{print()}, \code{fit()} and \code{predict()}-methods are avaiable.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'summary(fit)
#'@export
summary.mquantfit <- function(object, ...){
    cat("Call:\n")
    dput(object$call)
    cat("\n")

    # Creating a summary for each q-value estimated
    for (i.q in 1:length(object$q.values)) {
      cat("Estimation for m-quantile  q =",
          object$q.values[i.q],
          "with k = ",
          object$k,
          "\n")
      # Preparing and Printing Coefficient Matrix
      coef <- coefficients(object)
      if(object$method != "quasipoisson" & object$method != "nb") {
      out.q <- matrix(nrow = ncol(coef) - 1, ncol = 2)
      out.q[, 1] <- coef[i.q, 2:ncol(coef)]
      out.q[, 2] <- sqrt(object$var.beta[, i.q]) # Check if beta variance is implemented
      row.names(out.q) <- colnames(coef)[2:ncol(coef)]

      colnames(out.q) <- c("coefficients", "Std. Error")
      } else {
        out.q <- matrix(nrow = ncol(coef) - 1, ncol = 1)
        out.q[, 1] <- coef[i.q, 2:ncol(coef)]
        # out.q[, 2] <- sqrt(object$var.beta[, i.q]) # Check if beta variance is implemented
        row.names(out.q) <- colnames(coef)[2:ncol(coef)]

        colnames(out.q) <- c("coefficients")
      }
      cat("n =", length(object$Y), "\n")
      cat("\nCoefficients:\n")
      print(out.q)
      cat("Estimator for scale parameter sigma: ", object$scale[i.q], "\n")

      # Printing Residual Information
      cat("\nResiduals:\n")
      cat(t(summary(object$residuals[, i.q])), "\n")
      cat("\n")
      cat("Proportion of residuals smaller than 0:",
          mean(object$residuals[, i.q] < 0),
          "\n")
      cat("Residuals bounded at +-", object$k * object$scale[i.q], "\n")
      cat("Proportion of Huberised residuals:",
          mean(object$huberised.res[, i.q]) * 100,
          "%\n")

      if ((object$method == "continuous") &&
          (object$scale.estimator == "Mad"))
        cat("Pseudo R-squared (Bianchi et al, 2015): ",
            r2pseudo(object)[i.q],
            "\n")
      # Printing further Metrics
      cat("\n---\n")
        if(object$method == "continuous") {
          cat("Estimation with likelihood methodology\n")
          if(object$scale.estimator == "Mad")
            cat("Scale estimator: Mean Absolute Deviation (MAD)\n")
          if(object$scale.estimator == "cMad")
            cat("Scale estimator: corrected Mean Absolute Deviation (cMAD) by James Dawber (2017) \n")
          if(!((object$scale.estimator == "Mad") ||
             (object$scale.estimator == "cMad")))
            cat("Scale estimator: Custom scale estimator used\n")
        }


      if (object$method != "continuous"){
        cat("Estimation with Robust Quasi-Likelihood Estimation\n")
      }


      if (length(object$q.values) > 1 &&
          i.q != length(object$q.values)) {
        cat("\n------------------------\n\n")
      }

    }
    if (!is.null(attr(object, "na.message")))
      cat(attr(object, "na.message"), "\n")
    invisible(object)
}



#'Show Huberised Residuals
#'@rdname huberised
#'@export
huberised <- function(x) {
  UseMethod("huberised")
}
#' Huberised Residuals
#'
#' The huberised-method returns a flag indicating if an observation has been down-weighted by the alghorithm. .
#'
#'@name huberised
#'@param x a mquantfit object, produced by calling the mquantreg-function.
#'@return a dataframe-object with booleans indicating if the observation has been Huberised.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit  <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'huberised(fit)
#'
#'@export
huberised.mquantfit <- function(x) {
    x$huberised.res
  }

#'Compute Q-Scores
#'@rdname qscores
#'@export
qscores <- function(x, ...) {
  UseMethod("qscores")
}

#' Compute Q-Scores
#'
#' The qscores()-method returns the qscores compute via inverse mquantile estimation. .
#'
#'@name qscores
#'@param x a mquantfit object, produced by calling the mquantreg-function.
#'@return a dataframe-object with the qscores of all observations.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <-  mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'qscores(fit)
#'
#'@export
qscores.mquantfit <- function(x, qgrid = c(0.1,
                                           0.15,
                                           0.2,
                                           0.25,
                                           0.30,
                                           0.35,
                                           0.4,
                                           seq(from = 0.45, to = 0.55, by = 0.005),
                                           0.60,
                                           0.65,
                                           0.70,
                                           0.75,
                                           0.8,
                                           0.85,
                                           0.9)) {
  if (x$compute.qscores) {
      x$q.scores
    } else {
      x$call$compute.qscores <- TRUE
      x$call$qgrid <- qgrid
      a <- eval(x$call)
      a$q.scores
    }
  }


#' Plotting M-Quantile Diagnostic Plots
#'
#' The plot-method creates different diagnostic plots for m-quantile regression.
#'
#'@param x a mquantfit object, produced by calling the mquantreg-function.
#'@param plottype string, which specifies the plot to produce.
#'#'The avaiable plot types are: \cr
#'\itemize{
#'\item \code{"fit.res"} - Fitted vs. Residuals Plot
#'\item \code{"prop.res"} - Proportion of Residuals Smaller than 0 vs. Fitted Values
#'\item \code{"fitted.observed"} - Fitted vs. Observerd Values
#'\item \code{"optk"} - Comparing the Residual Distribution with Normal Distribution
#'\item \code{"coef.q"} - Coefficients over q
#'\item \code{"prop.huber.k"} - Proportion of huberised residuals over k
#'\item \code{"prop.huber.q"} - Proportion of huberised residuals over q
#'}
#'@param add.qvals numeric, containing additional m-quantiles to be shown in plot.
#'@param show.huberised boolean controlling if the plot should show the squared residuals instead.
#'@param sqr.res boolean controlling if the plot should show the squared residuals instead.
#'@param  kstart numeric, the first k estimated and shown in the plot.
#'@param  kend  numeric, the last k estimated and shown in the plot.
#'@param  kstep numeric, controls the steps in which k is increased between kstart and kend.
#'@param  coef numeric, choose index of coefficient to plot.
#'@param ...	Arguments to be passed to methods, such as graphical parameters.
#'@return a dataframe-object with the qscores of all observations.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "prop.huber.k", kstart = 0.1, kend = 4)
#'
#'@export
plot.mquantfit <- function(x,
                           plottype = "fit.res",
                           add.qvals = FALSE,
                           show.huberised = TRUE,
                           sqr.res = FALSE,
                           kstart = 0.5,
                           kend = 3,
                           kstep = 0.1,
                           coef = FALSE,
                           ...) {
  ##############################################
  ## Check if valid method is selected
  if (!(plottype %in% c("fit.res",
                        "prop.res",
                        "fitted.observed",
                        "optk",
                        "coef.q",
                        "prop.huber.k",
                        "prop.huber.q"
                        )
        )) {
    stop("Chosen method is not implemented. Avaiable methods are: fit.res, prop.res,
         fitted.observed, optk, coef.q, prop.huber.k, prop.huber.q"
         )
  }

  switch(plottype,
         fit.res = ggplot.fit.res(x, sqr.res = sqr.res, add.qvals = add.qvals),
         prop.res = ggplot.prop.res(x, add.qvals = add.qvals),
         fitted.observed = ggplot.fitted.observed(x, add.qvals = add.qvals),
         optk = ggplot.optk(x,
                            kstart = kstart,
                            kend = kend,
                            kstep = kstep
                            ),
         coef.q = ggplot.coef.q(x, coef),
         prop.huber.k = ggplot.prop.huber.k(x,
                                            kstart = kstart,
                                            kend = kend,
                                            kstep = kstep,
                                            add.qvals = add.qvals
                                            ),
         prop.huber.q = ggplot.prop.huber.q(x)
         )
}

#' Pseudo R-Squared for Continuous M-Quantile Estimation
#'@rdname r2pseudo
#'@export
r2pseudo <- function(fit) {
  UseMethod("r2pseudo")
}

#' Pseudo R-Squared for Continuous M-Quantile Estimation
#'
#' Computes the pseudo R-Squared for continuous M-Quantile Models with MAD scale estimator. The procedure was proposed by Bianchi (2015).
#' Bianchi, A., Fabrizi, E., Salvati, N. and Tzavidis, N. (eds) (2015). M-quantile re-
#  gression: diagnostics and parametric representation of the model.
#'@name r2pseudo
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@return numeric, the pseudo r-square value.
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'r2pseudo(fit)
#'
#'@export
r2pseudo.mquantfit <- function(fit) {
  ##############################################
  ## Check if model is continuous and the MAD scale estimator is used
  if ((fit$scale.estimator == "Mad") &&
      (fit$method == "continuous") && (fit$k > 0)) {
    q     <-  fit$q.values
    k     <-  fit$k
    scale <-  fit$scale
    k     <- fit$k
    r2.q  <- NULL
    # Compute pseudo R-Square
    # Bianchi, A., Fabrizi, E., Salvati, N. and Tzavidis, N. (eds) (2015). M-quantile re-
    # gression: diagnostics and parametric representation of the model.
    huber.loss <- function(x, k = 1.345) {
      ifelse (abs(x) <= k, 0.5 * x ^ 2, k * (abs(x) - 0.5 * k))
    }

    #Computing restricted Model
    for(i in 1:length(q)){
      q_i <- q[i]
      # calculate full and reduced model
      tmp_full <- mquantreg(formula = fit$formula, data = fit$data, q = q_i,
                            method = "continuous", k = k)
      tmp_full$call$data    <- data.frame("y" = fit$Y)
      tmp_full$call$formula <- "y ~ 1"
      tmp_res  <- eval(tmp_full$call)

      num <- sum(apply(tmp_full$residuals / tmp_full$scale, 1, huber.loss, k = k))
      den <- sum(apply(tmp_res$residuals / tmp_res$scale,1,
                       huber.loss,
                       k = k))
      r2.q[i] <- (1 - (num / den))
    }
    # fit$call$data    <- data.frame("y"  = fit$Y)
    # fit$call$formula <- "y ~ 1"
    # fit_restricted   <- eval(fit$call)
    # return(1 - (den / num))
    return(r2.q)
  } else if ((fit$scale.estimator == "Mad") &&
             (fit$method == "continuous") && (fit$k > 0)){
    p <- ncol(fit$X)
    q <- fit$q.values
    n <- nrow(fit$residuals)
    q <- fit$q.values
    k <- fit$k
    r2.q <- NULL


    for (i in 1:length(q)) {
      #i <- 1
      q_i <- q[i]
      # calculate full and reduced model
      tmp_full <- mquantreg(formula = fit$formula, data = fit$data, q = q_i,
                            method = "continuous", k = k)

      dat_tmp    <- data.frame("y" = tmp_full$Y)
      tmp_red    <- mquantreg(formula = y ~ 1, data = dat_tmp, q = q_i,
                              method = "continuous", k = k)

      res_full <- tmp_full$resid
      res_red  <- tmp_red$resid

      sm              <- tmp_full$scale
      w               <- (k*(abs(res_full)-k/2))*((res_full/sm)<=-k | (res_full/sm)>= k)+0.5*res_full^2*(-k<(res_full/sm) & (res_full/sm)<k )
      ww              <- 2 * (1 - q[i]) * w
      ww[res_full> 0] <- 2 * q[i] * w[res_full > 0]
      w               <- ww
      V_full          <- sum(w)

      sm2             <- sm
      w2              <- (k*(abs(res_red)-k/2))*((res_red/sm2)<=-k | (res_red/sm2)>= k)+0.5*res_red^2*(-k<(res_red/sm2) & (res_red/sm2)<k )
      ww2             <- 2 * (1 - q[i]) * w2
      ww2[res_red> 0] <- 2 * q[i] * w2[res_red > 0]
      w2              <- ww2
      V_red           <- sum(w2)

      r2.q[i]         <- 1 - (V_full / V_red)

    }

    return(r2.q)

  } else {
    stop("The pseudo R-squared can only be computed for continuous data using the MAD scale estimator.")
    den = NULL
    num = NULL
  }

}

#' M-Quantile Estimation for Continuous Data using MLE.
#'
#' The function implements m-quantile estimation for continuous data first proposed by Breckling and Chambers (1988).
#' This implementation uses the paremetric representation introduced by Bianchi et al. (2015).
#' The function is called by the mquantreg-function, when the method "continuous" is selected. Important:
#' Standalone use is possible, but not advised.
#'
#'@param x matrix, containing independent variables.
#'@param y vector, containing the continuous dependent variable.
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'@param k a number greater than 0. k is the parameter for the huber psi-function, the default value is 1.345.
#'@param scale.estimator string, defines which function is used to estimate the scale parameter.
#'Following function are implemented in the package: \cr
#'\itemize{
#'\item \code{"Mad"} - Mean Absolute Deviation (default)
#'\item \code{"cMad"} - Corrected Mean Absolute Deviation
#'}
#'@param maxit an integer, which defines the maximal number of iteration
#'before the algorithm stops although it has not converged yet.
#'If the maximal number of iteration is reached a warning message will be prompted.
#'@param acc defines convergence criteria
#'@param case.weights vector of observation weights;
#'if supplied, the algorithm fits to minimize the sum of the weights multiplied into the absolute residuals.
#'The length of weights must be the same as the number of observations.
#'The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous
#'@param var.weights dataframe, if supplied the residuals are scaled.
#'@param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting.
#'This should be NULL or a numeric vector of length equal to the number of cases.
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@param compute.qscores boolean, controls if q-scores are estimated.
#'@return See \code{\link{mquantreg}} for details.
#'\code{summary()}, \code{print()}, \code{fitted()} and \code{predict()}-methods are avaiable.
#'@references
#'\itemize{
#'\item Breckling, J., & Chambers, R. (1988). M-quantiles. Biometrika, 75(4), 761-771.
#'\item Bianchi, A., Fabrizi, E., Salvati, N., & Tzavidis, N. (2015).
#'M-quantile Regression: Diagnostics and the parametric representation of the model.
#'Book of Abstracts SIS - CLADAG, 303-306
#'\item Bianchi, A., & Salvati, N. (2015). Asymptotic properties and variance estimators
#'of the M-quantile regression coefficients estimators. Commun.Stat. Theory., 44, 2016–2429.
#'}
#'
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                    real.betas = c(0.1, 0.3, 0.1 ),
#'                    response.type = "continuous.normal",
#'                    measurement.error = 0.05)
#'
#'fit <- mquantreg(formula = Y ~ x1 + x2,data = df, q  = 0.5, method = "continuous")
#'print(fit)
#'
mqcont <- function(x,
                   y,
                   q = 0.5,
                   k = 1.345,
                   scale.estimator = "Mad",
                   offset = NULL, # rep(0, nrow(x))
                   case.weights = NULL, # rep(1, nrow(x))
                   var.weights = NULL, # rep(1, nrow(x))
                   maxit = 1000,
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
                             0.9)){

  if (is.null(offset)) {
    offset <- rep(0, nrow(x))
  }

  if (is.null(case.weights)){
    case.weights <- rep(1, nrow(x))
  }


  if (is.null(var.weights)) {
    var.weights <- rep(1, nrow(x))
  }

  #probability density function M-quantile

  if(scale.estimator == "cMad"){
    scale.estimator.fnc <- function(resid, var.weights ) {
      median(abs(resid / sqrt(var.weights) - median(resid / sqrt(var.weights)))) / 0.6745 # qnorm(3 / 4)
    }
  }

  if(scale.estimator == "Mad"){
    scale.estimator.fnc <- function(resid, var.weights) {
      median(abs(resid/sqrt(var.weights))) / 0.6745 # qnorm(3 / 4)
    }
  }

  pdf_mqY <- function(y, m, s, q, k) {

    myrho_huber <- function(y = y, m = m, q = q, k = k, s = s){
      e         <- (y - m) / s
      w         <- ((e^2 / 2) * (-k <= e & e <= k) + (k * abs(e) - k^2 / 2) *
                      (e < -k | e > k))
      ww        <- 2 * (1 - q) * w
      ww[e > 0] <- 2 * q * w[e > 0]
      w <- ww
      w
    }

    Bq <- function(q = q, k = k, s = s) {
      s * (sqrt(pi / q) * (pnorm(k * sqrt(2 * q)) - 0.5) + sqrt(pi / (1 - q)) *
             (pnorm(k * sqrt(2 * (1 - q))) - 0.5) + (1 / (2 * k * q)) * exp(-q * k^2) +
             (1 / (2 * k * (1 - q))) * exp(-(1 - q) * k^2))
    }

    pdf <- (1 / Bq(q = q , k = k, s = s)) * exp( -myrho_huber(y = y, m = m,q = q ,k = k, s = s))
    pdf
  }

  LLL <- function(par) {
    s     <- par[1]
    n     <- length(y)
    phiik <- pdf_mqY(y = y, m = mik, s = s, q = q[i], k = k)
    fik   <- array(ifelse(phiik > 0, phiik, 1e-005), c(n, 1))

    # calculation of log-likelihood
    lik    <- fik
    li     <- log(lik)
    likeli <- sum(li)

    -likeli
  }


  cs <- function(init, step = 10, maxit = 500, tol = 1e-5,...) { # tol = 1e-5
    b <- init
    e <- cbind(diag(length(b)), -diag(length(b)))
    for(i in 1:maxit) {
      flag = FALSE
      for(j in 1:ncol(e)) {
        if(LLL(par = b + e[,j] * step,...) < LLL(par=b, ...)) {
          b <- b + e[,j] * step
          flag <- TRUE
          break
        }
      }
      if(!flag) {
        step <- step/2
        if(step < tol) break
      }
    }
    if(i == maxit){
      print(paste("Convergence was not achieved in", maxit, "iterations."))
    }
    return(b)
  }



  n <- length(y)
  p <- ncol(as.matrix(x))

  qest        <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt         <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit        <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres        <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar        <- matrix(0, nrow = ncol(x), ncol = length(q))
  huberised   <- matrix(1, nrow = nrow(x), ncol = length(q))

  qvar.matrix <- array(rep(0, ncol(x) * ncol(x)), dim = c(ncol(x), ncol(x), length(q)))
  qscale      <- NULL
  qloglike    <- NULL
  qAIC        <- NULL

  for(i in 1:length(q)) {

    init  <- QRLM(x = x, y = y, k = k, maxit = maxit, acc = acc, q = q[i],
                  scale.estimator = scale.estimator)

    resid <- y - x %*% init$coef
    diff  <- 1
    niter <- 1
    beta  <- c(init$coef)
    scale <- init$qscale

    while(diff > acc) {
      scale.old <- scale

      for (iiter in 1:maxit) {
        w             <- psi.huber(resid / (scale.old * sqrt(var.weights)), k = k) *
                          case.weights
        huberised[, i] <- (psi.huber((resid) / scale * sqrt(var.weights), k = k) == 1)
        ww            <- 2 * (1 - q[i]) * w
        ww[resid > 0] <- 2 * q[i] * w[resid > 0]
        w             <- ww
        temp          <- lm.wfit(x, y - offset, w, method = "qr")
        resid         <- temp$residuals
        convi         <- sum(abs((temp$coef - beta) / (beta)))
        beta          <- temp$coef

        done          <- (convi <= acc)
        if (done)
          break
      }

      qest[,i] <- temp$coef
      beta     <- temp$coef
      mik      <- offset + temp$fitted
      c.optim  <- cs(init = c(scale.old), step = 0.1)
      scale    <- c.optim[1]
      diff     <- abs((scale - scale.old) / (scale.old))
      niter    <- niter + 1
      done1    <- (niter > maxit)
      if (done1)
        break
    }

    qscale[i]  <- scale

    if (done1)
      warning(paste("mquantreg failed to converge in", maxit, "steps at q = ", q[i]))

    resid    <- temp$residuals
    qres[,i] <- resid
    qfit[,i] <- offset + temp$fit
    qwt[,i]  <- w

    # estimation of variance for regression coefficients
    resid1   <- resid / scale
    Epsi2    <- (sum((w * resid1)^2)/(n - p))
    Epsi     <- (1 / scale) * (sum(2 * (q[i] * (0 <= resid1 & resid1 <= k) + (1 - q[i]) *
                (-k <= resid1 & resid1 < 0))) / n)
    xt       <- t(x)
    var.beta <- (((Epsi2) / Epsi^2) * solve(xt %*% x))
    qvar[,i] <- diag(var.beta)
    qvar.matrix[,,i] <- var.beta

    # AIC
    qloglike[i] <- LLL(c(scale))
    qAIC[i]     <- 2 * (p + 2) + 2 * qloglike[i]


    sqrt.beta <- NULL
    for (i in 1:p) {
      sqrt.beta[i] <- sqrt(var.beta[i, i])
    }

    q.scores <- NULL
    # compute q-values
    if((compute.qscores) && k > 0.01 ) {
      q.scores <- mqcont.qscores(x               = x,
                                 y               = y,
                                 offset          = offset,
                                 case.weights    = case.weights,
                                 var.weights     = var.weights,
                                 maxit           = maxit,
                                 acc             = acc,
                                 qgrid           = qgrid,
                                 k.value         = k,
                                 epsilon         = 0.01,
                                 scale.estimator = scale.estimator
                                 )$qscores
    }


  }


  list(
    fitted.values   = qfit,
    residuals       = qres,
    q.values        = q,
    q.weights       = qwt,
    coefficients    = qest,
    scale           = qscale,
    q.scores        = q.scores,
    huberised.res   = !(huberised),
    var.beta        = qvar,
    var.beta.matrix = qvar.matrix
  )

}

QRLM <- function (x,
                  y,
                  case.weights = rep(1, nrow(x)),
                  var.weights = rep(1, nrow(x)),
                  scale.estimator,
                  offset = rep(0, nrow(x)),
                  k = 1.345,
                  maxit = 20,
                  acc = 1e-04,
                  q = 0.5) {

  irls.delta <- function(old, new) {
    sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  }

  if(scale.estimator == "cMad"){
    scale.estimator.fnc <- function(resid, var.weights ) {
      median(abs(resid / sqrt(var.weights) - median(resid / sqrt(var.weights)))) / 0.6745
    }
  }

  if(scale.estimator == "Mad"){
    scale.estimator.fnc <- function(resid, var.weights) {
      median(abs(resid/sqrt(var.weights))) / 0.6745
    }
  }

  # Compute initial weights
  init.weights <- rep(1, nrow(x))
  w            <- (init.weights * case.weights) / var.weights

  # compute initial estimators
  init.fit   <- lm.wfit(x, y - offset, w, method = "qr")
  coef       <- init.fit$coef
  resid      <- init.fit$resid

  resid.init <- resid
  done       <- FALSE
  conv       <- NULL

  scale      <- scale.estimator.fnc(resid, var.weights)
  qest       <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt        <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit       <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres       <- matrix(0, nrow = nrow(x), ncol = length(q))
  huberised  <- matrix(1, nrow = nrow(x), ncol = length(q))
  qscale     <- numeric(length(q))


  for(i in seq_along(q)) {
    for (iiter in 1:maxit) {
        resid.old <- resid
        scale     <- scale.estimator.fnc(resid, var.weights)

        if (scale == 0) {
          done <- TRUE
          break
        }

      psi.res        <- psi.huber(resid / (scale * sqrt(var.weights)), k = k)
      huberised[, i] <- (psi.huber((resid) / scale * sqrt(var.weights), k = k) == 1)
      psi.res        <- (psi.res * case.weights) / var.weights
      w              <- 2 * (1 - q[i]) * psi.res
      w[resid > 0]   <- 2 * q[i] * psi.res[resid > 0]

      iter.fit <- lm.wfit(x, y - offset, w, method = "qr")

      coef  <- iter.fit$coef
      resid <- iter.fit$residuals
      convi <- irls.delta(resid.old, resid)
      conv  <- c(conv, convi)
      done  <- (convi <= acc)

      if (done)
        break
    }
    if (!done)
      warning(paste("QRLM failed to converge in", maxit, "steps at q = ", q[i]))

    qest[, i]  <- coef
    qscale[i]  <- scale
    qwt[, i]   <- w
    qfit[, i]  <- offset + iter.fit$fitted.values
    qres[, i]  <- resid
  }
  list(fitted.values = qfit,
       residuals     = qres,
       q.values      = q,
       q.weights     = qwt,
       coefficients  = qest,
       qscale        = qscale,
       huberised.res = !(huberised)
  )
}


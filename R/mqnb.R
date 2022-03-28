#' M-Quantile Estimation for Count Data with Negative Binomial Distribution
#'
#' The function implements m-quantile estimation for count data as proposed by Chambers et al.(2014).
#' The algorithm is based on robust quasi-likelihood estimation by Cantoni & Ronchetti (2001).
#' The function is called by the mquantreg-function,
#' when the method "nb" is selected. Important: Standalone use is possible, but not advised.
#'
#'@param x matrix, containing independent variables
#'@param y matrix, containing dependent variables
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'@param k a number greater than 0. k is the parameter for the huber psi-function, the default value is 1.6 for standalone use.
#'@param maxit an integer, which defines the maximal number of iteration
#'before the algorithm stops allthough it has not converged yet.
#'If the maximal number of iteration is reached a warning message will be prompted.
#'@param acc defines convergence criteria
#'@param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting.
#'This should be NULL or a numeric vector of length equal to the number of cases.
#'@param case.weights vector of observation weights;
#'if supplied, the algorithm fits to minimize the sum of the weights multiplied into the absolute residuals.
#'The length of weights must be the same as the number of observations.
#'The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous.
#'@param weights.x vector of weights, for estimation of robust glm to obtain initial weights in the discrete algorithms.
#'@param var.weights dataframe, if supplied the residuals are scaled.
#'@param theta.init numeric, initial theta.
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@param compute.qscores boolean, controls if q-scores are estimated.
#'@return See \code{\link{mquantreg}} for details.
#'\code{summary()}, \code{print()}, \code{fit()} and \code{predict()}-methods are avaiable.
#'@references
#'Cantoni, E., & Ronchetti, E. (2001). Robust inference for generalized linear models. Journal of the American Statistical Association, 96(455), 1022-1030.
#'Chambers, R., Dreassi, E., & Salvati, N. (2014). Disease mapping via negative binomial regression M‐quantiles. Statistics in medicine, 33(27), 4805-4824.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "count.nb",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = Y ~ x1 + x2, data = df, q  = 0.5, method = "nb")
#'print(fit)

mqnb <- function(x,
                 y,
                 q = 0.5,
                 k = 1.6,
                 maxit = 1000,
                 acc = 1e-04,
                 offset = NULL, # rep(1, nrow(x))
                 case.weights = NULL, # rep(1, nrow(x))
                 weights.x = FALSE,
                 var.weights = NULL, # rep(1, nrow(x))
                 theta.init = 1,
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
                          0.9),
                 compute.qscores = FALSE
                 ) {

  if (is.null(offset)) {
    offset <- rep(1, nrow(x))
  }

  if (is.null(case.weights)){
    case.weights <- rep(1, nrow(x))
  }


  if (is.null(var.weights)) {
    var.weights <- rep(1, nrow(x))
  }

  #estimation of theta
  theta.mm.mq = function(y,
                         fit,
                         q,
                         theta,
                         dfr,
                         k = 1.345,
                         acc = 0.001,
                         maxit = 20) {
    ExpectedZ <- function(y, mu, q, theta, k) {
      V       <- mu + (mu ^ 2) / theta
      r.stand <- (y - mu) / sqrt(V)

      #we compute i1 and i2
      jinf <- floor(mu - k * sqrt(V))
      jsup <- floor(mu + k * sqrt(V))
      w.q  <- 2 * (q * (r.stand > 0) + (1 - q) * (r.stand <= 0))
      Ez   <- (k ^ 2) * ((pnbinom(jsup, size = theta, mu = mu, lower.tail = F) +
                pnbinom((jinf), size = theta, mu = mu))) + (mu / V) *
                (dnbinom(jinf, size = theta, mu = mu) * (jinf / theta) *
                (theta + 1 + jinf) - dnbinom(jsup, size = theta, mu = mu) *
                (jsup / theta) * (theta + 1 + jsup) +
                (pnbinom((jsup - 1), size = theta, mu = mu) -
                pnbinom(jinf, size = theta, mu = mu))) + (mu ^ 2 / V) *
                (dnbinom(jinf, size = theta, mu = mu) *
                ((jinf - jinf * theta - theta ^ 2) / theta ^ 2) -
                dnbinom(jsup, size = theta, mu = mu) *
                ((jsup - jsup * theta - theta ^ 2) / theta ^ 2) +
                (1 / theta) * (pnbinom((jsup - 1), size = theta, mu = mu) -
                pnbinom(jinf, size = theta, mu = mu)))
      mean(w.q ^ 2 * Ez)
    }

    basepsiq <- function(res, q) {
      hub.res.sym      <- psi.huber(res, k = k) * res
      hub.res          <- 2 * (1 - q) * hub.res.sym
      hub.res[res > 0] <- 2 * q * hub.res.sym[res > 0]
      return(hub.res)
    }


    diff      <- 1
    theta.old <- theta
    resid     <- (y - fit)
    iter.cnt  <- 0

    while (diff > acc) {
      iter.cnt <- iter.cnt + 1
      V        <- fit + (fit ^ 2) / theta.old
      res.std  <- resid / (sqrt(V))
      EZ       <- ExpectedZ(y, mu = fit, q = q, theta = theta.old, k = k)
      hub.res  <- basepsiq(res = res.std, q = q)
      wi0      <- hub.res / res.std

      theta.new <- (1 / dfr) * sum((wi0 ^ 2 * ((y - fit) / sqrt(theta.old *
                    fit + fit ^ 2)) ^ 2) / EZ)
      theta.new <- 1 / theta.new
      diff      <- (theta.new - theta.old) ^ 2
      theta.old <- theta.new
      crit1     <- (theta.new > 100000)
      crit2     <- (iter.cnt > maxit)
      if (crit1 | crit2) {
        break
      }

    }
    theta.new
  }

  #Stopping rule
  irls.delta <- function(old, new) {
    abs(max(old - new)) / abs(max(old))
  }

  n <- length(case.weights)
  if (weights.x) {
    w.x <- sqrt(1 - hat(x))
  } else {
    w.x <- rep(1, length = n)
  }

  #We fit the glm.rob for computing the starting values
  p <- ncol(x)

  done      <- FALSE
  conv      <- NULL
  qest      <- matrix(0, nrow = ncol(x), ncol = length(q))
  qfit      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar      <- matrix(0, nrow = ncol(x), ncol = length(q))
  huberised <- matrix(1, nrow = nrow(x), ncol = length(q))
  qtheta    <- NULL

  for (i in 1:length(q)) {
    #We define the starting values
    temp.mq.poisson <- mqpoisson(x            = x,
                                 y            = y,
                                 offset       = offset,
                                 case.weights = case.weights, # rep(1, nrow(x))
                                 maxit        = maxit,
                                 acc          = acc,
                                 weights.x    = weights.x,
                                 q            = q[i],
                                 k            = k
                                 )

    theta <- theta.mm.mq(y,
                         fit   = temp.mq.poisson$fitted.values,
                         q     = q[i],
                         theta = theta.init,
                         dfr   = (n - p),
                         k     = k,
                         acc = acc,
                         maxit = maxit
                         )

    resid <- y - temp.mq.poisson$fitted.values
    fit   <- temp.mq.poisson$fitted.values
    a.j   <- case.weights
    w     <- case.weights
    coef  <- temp.mq.poisson$coef

    for (iiter in 1:maxit) {
      resid.old <- resid
      coef.old  <- coef

      # We define the probability mu=t*exp(xb)
      probab   <- fit
      mu       <- probab
      deriv.mu <- mu

      #We define the variance
      V <- probab + (probab ^ 2) / theta

      #We define the scale
      scale <- c(sqrt(V))

      #We standardize the residuals
      r.stand <- (y - mu) / sqrt(V)

      #we compute i1 and i2
      jinf <- floor(mu - k * sqrt(V))
      jsup <- floor(mu + k * sqrt(V))

      #We compute the values of a_j(b)
      if (k == Inf) {
        a.j <- rep(1, n)
      }

      if (k != Inf) {
        a.j <- (-k) * pnbinom(jinf, size = theta, mu = mu) + k *
                (1 - pnbinom((jsup), size = theta, mu = mu)) + mu / sqrt(V) *
                (pnbinom(jinf, size = theta, mu = mu) -
                pnbinom((jinf - 1), size = theta, mu = mu)) * (1 + jinf / theta) -
                mu / sqrt(V) * (pnbinom(jsup, size = theta, mu = mu) -
                pnbinom((jsup - 1), size = theta, mu = mu)) * (1 + jsup / theta)
      }

      a.j <- 2 * a.j * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))


      #we define a part of w_j
      w <- diag(c(mu) / scale) * diag(c(w.x))

      #we compute psi_q(res)
      hub.res.sym        <- psi.huber((resid) / scale, k = k) * case.weights *
                            ((resid) / scale)
      huberised[, i]     <- (psi.huber((resid) / scale, k = k) == 1)
      hub.res            <- 2 * (1 - q[i]) * hub.res.sym
      hub.res[resid > 0] <- 2 * q[i] * hub.res.sym[resid > 0]

      #we compute psi_q(r )-E(psi_q(r ))
      A <- (hub.res - a.j)


      if (k == Inf) {
        esp.carre.cond <- rep(1, n)
      }

      if (k != Inf) {
        esp.carre.cond <- k * (mu / V) * ((pnbinom(jinf, size = theta, mu = mu) -
                           pnbinom((jinf - 1), size = theta, mu = mu)) *
                           ((theta + jinf) / theta) +
                           (pnbinom(jsup, size = theta, mu = mu) -
                           pnbinom((jsup - 1), size = theta, mu = mu)) *
                           ((theta + jsup) / theta)) + (mu / V ^ (3 / 2)) *
                           ((pnbinom(jinf, size = theta, mu = mu) -
                           pnbinom((jinf - 1), size = theta, mu = mu)) *
                           (jinf / theta) * (theta + 1 + jinf) -
                           (pnbinom(jsup, size = theta, mu = mu) -
                           pnbinom((jsup - 1), size = theta, mu = mu)) *
                           (jsup / theta) * (jsup + 1 + theta) +
                           pnbinom((jsup - 1), size = theta, mu = mu) -
                           pnbinom(jinf, size = theta, mu = mu)) +
                           (mu ^ 2 / V ^ (3 / 2)) * ((pnbinom(jinf, size = theta, mu = mu) -
                           pnbinom((jinf - 1), size = theta, mu = mu)) *
                           ((jinf - jinf * theta - theta ^ 2) / theta ^ 2) -
                           (pnbinom(jsup, size = theta, mu = mu) -
                           pnbinom((jsup - 1), size = theta, mu = mu)) *
                           ((jsup - jsup * theta - theta ^ 2) / theta ^ 2) +
                           (1 / theta) * (pnbinom((jsup - 1), size = theta, mu = mu) -
                           pnbinom(jinf, size = theta, mu = mu)))
      }

      b.j <- 2 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
      B   <- diag(c(mu * b.j))

      #We estimate betas
      temp  <- coef + solve(t(x) %*% w %*% B %*% x) %*% t(x) %*% w %*% A
      coef  <- temp
      eta   <- x %*% coef
      fit   <- offset * exp(eta)
      resid <- y - fit

      convi <- irls.delta(coef.old, coef)
      conv  <- c(conv, convi)
      done  <- (convi <= acc)

      if (done) {
        break
      }
    }

    if (!done) {
      warning(paste("mqnb failed to converge in", maxit, "steps at q = ", q[i]))
    }

    # Asymptotic estimated variance of the robust estimator

    probab   <- fit
    mu       <- probab
    deriv.mu <- mu

    #We define the variance
    V <- probab + (probab ^ 2) / theta

    #We define the scale
    scale <- c(sqrt(V))

    #We standardize the residuals
    r.stand <- (y - mu) / sqrt(V)


    if (k == Inf) {
      esp.cond <- rep(1, n)
    } else {
      esp.cond <- (-k) * pnbinom(jinf, size = theta, mu = mu) + k *
                (1 - pnbinom((jsup), size = theta, mu = mu)) + mu / sqrt(V) *
                (pnbinom(jinf, size = theta, mu = mu) -
                pnbinom((jinf - 1), size = theta, mu = mu)) * (1 + jinf / theta) -
                mu / sqrt(V) * (pnbinom(jsup, size = theta, mu = mu) -
                pnbinom((jsup - 1), size = theta, mu = mu)) * (1 + jsup / theta)
    }

    esp.cond <- 2 * esp.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
    a.const  <- apply(x * as.vector(1 / n / sqrt(V) * w.x * esp.cond * deriv.mu), 2, sum)

    if (k == Inf) {
      esp.carre.cond <- rep(1, n)
    } else {
      esp.carre.cond <- (k ^ 2) * (1 - (pnbinom(jsup, size = theta, mu = mu) - pnbinom(jinf, size = theta, mu = mu))) +
                        mu / V * ((pnbinom(jinf, size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)) *
                        jinf / theta * (theta + 1 + jinf) - (pnbinom(jsup, size = theta, mu = mu) - pnbinom((jsup - 1), size = theta, mu = mu) *
                        jsup / theta * (theta + 1 + jsup)) + (pnbinom((jsup - 1), size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu))) +
                        (mu ^ 2) / V * ((pnbinom(jinf, size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)) *
                        ((jinf - jinf * theta - theta ^ 2) / theta ^ 2) - (pnbinom(jsup, size = theta, mu = mu) - pnbinom((jsup - 1), size = theta, mu = mu) *
                        ((jsup - jsup * theta - theta ^ 2) / theta ^ 2)) + 1 / theta * (pnbinom((jsup - 1), size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)))
    }

    esp.carre.cond <- 4 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0)) ^ 2
    matQaux <- as.vector(esp.carre.cond / V * w.x ^ 2 * deriv.mu ^ 2)
    matQ1   <- (1 / n) * t(x) %*% (matQaux * x)
    matQ2   <- a.const %*% t(a.const)
    matQ    <- matQ1 - matQ2

    if (k == Inf) {
      esp.psi.score = 1 / sqrt(V)
    } else {
      esp.psi.score <- mu * k / V * ((pnbinom(jinf, size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)) *
                       (jinf + theta) / theta + (pnbinom(jsup, size = theta, mu = mu) - pnbinom((jsup - 1), size = theta, mu = mu)) *
                       (jsup + theta) / theta) + mu / V ^ (3 / 2) * ((pnbinom(jinf, size = theta, mu = mu) -
                       pnbinom((jinf - 1), size = theta, mu = mu)) * jinf / theta * (theta + 1 + jinf) -
                       (pnbinom(jsup, size = theta, mu = mu) - pnbinom((jsup - 1), size = theta, mu = mu) *
                       jsup / theta * (theta + 1 + jsup)) + (pnbinom((jsup - 1), size = theta, mu = mu) -
                       pnbinom((jinf - 1), size = theta, mu = mu))) + (mu ^ 2) / V ^ (3 / 2) *
                       ((pnbinom(jinf, size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)) *
                       ((jinf - jinf * theta - theta ^ 2) / theta ^ 2) - (pnbinom(jsup, size = theta, mu = mu) -
                       pnbinom((jsup - 1), size = theta, mu = mu) * ((jsup - jsup * theta - theta ^ 2) / theta ^ 2)) +
                       1 / theta * (pnbinom((jsup - 1), size = theta, mu = mu) - pnbinom((jinf - 1), size = theta, mu = mu)))
    }

    esp.psi.score <- 2 * esp.psi.score  * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
    matMaux <- as.vector(esp.psi.score / sqrt(V) * w.x * deriv.mu ^ 2)
    matM    <- 1 / n * t(x) %*% (matMaux * x)
    matMinv <- solve(matM)

    as.var  <- 1 / n * matMinv %*% matQ %*% matMinv

    qest[, i] <- coef
    qfit[, i] <- fit
    qres[, i] <- y - fit
    qtheta[i] <- theta
    qvar[, i] <- as.numeric(round(diag(as.var), 4))
  }

  if (compute.qscores) {
    q.scores <- mqnb.qscores(x,
                             y,
                             offset = offset,
                             case.weights = case.weights,
                             maxit = maxit,
                             acc = acc,
                             weights.x = weights.x,
                             qgrid = qgrid,
                             k = k,
                             epsilon = 0.01)$qscores
    } else {
    q.scores = NULL
    }

  list(fitted.values = qfit,
       var.beta      = qvar,
       residuals     = qres,
       q.values      = q,
       coefficients  = qest,
       q.theta       = qtheta,
       scale         = scale,
       q.scores      = q.scores,
       huberised.res = !(huberised)
       )
}

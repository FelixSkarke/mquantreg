#' M-Quantile Estimation for Count Data with Poisson Distribution
#'
#' The function implements m-quantile estimation for count data as proposed by Tzavidis et al. (2015).
#' The algorithm is based on robust quasi-likelihood estimation by Cantoni & Ronchetti (2001).
#' The function is called by the mquantreg-function, when method "poisson" is selected.
#' Important: Standalone use is possible, but not advised.
#'
#'@param x matrix, containing independent variables.
#'@param y vector, containing dependent variable.
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'Can also be a vector of m-quantiles.
#'@param k a number greater than 0. k is the parameter for the huber psi-function and
#'defaults to 1.6 for standalone use.
#'@param maxit a integer, which defines the maximum number of iterations
#'before the estimation process is halted. If the maximum number of iterations is
#'reached a warning message will be prompted.
#'If the maximal number of iteration is reached a warning message will be prompted.
#'@param acc defines convergence criteria.
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
#'@param var.weights dataframe, if supplied the residuals are scaled.
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@param compute.qscores boolean, controls if qscores are to be computed.Defaults to
#'\code{FALSE}.
#'@return See \code{\link{mquantreg}} for details.
#'\code{summary()}, \code{print()}, \code{fitted()} and \code{predict()}-methods are avaiable.
#'@references
#'\itemize{
#'\item Cantoni, E., & Ronchetti, E. (2001). Robust inference for generalized linear models. Journal of the American Statistical Association, 96(455), 1022-1030.
#'\item Tzavidis, N., Ranalli, M. G., Salvati, N., Dreassi, E., & Chambers, R. (2015). Robust small area prediction for counts. Statistical methods in medical research, 24(3), 373-395.
#'}
#'@examples
#'
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "count.poisson",
#'                    measurement.error = 0.05)
#'
#'fit <- mquantreg(formula = Y ~ x1 + x2,data = df, q  = 0.5, method = "poisson")
#'print(fit)
mqpoisson <- function(x,
                      y,
                      q = 0.5,
                      k = 1.6,
                      maxit,
                      acc = 1e-04,
                      offset = NULL,
                      case.weights = NULL,
                      weights.x = FALSE,
                      var.weights = NULL,
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
                      compute.qscores = FALSE) {

  if (is.null(offset)) {
    offset <- rep(1, nrow(x))
  }

  if (is.null(case.weights)){
    case.weights <- rep(1, nrow(x))
  }


  if (is.null(var.weights)) {
    var.weights <- rep(1, nrow(x))
  }

  #This is the robust glm with poisson for count data by Cantoni & Ronchetti
  glm.rob.poisson <- function(X,
                              y,
                              weights.x,
                              k = 1.6,
                              offset = offset) {

    mytol <- .Machine$double.eps ^ .25
    nb    <- length(y)

    basepsi <- function(x) {
      x * pmin(1, k / abs(x))
    }

    basepsiprime <- function(x) {
      1 * (abs(x) < k)
    }

    # Initializing....
    beta.old <- as.vector(glm(y ~ X - 1, family = poisson, offset = log(offset))$coeff)

    Xwith    <- X

    eta      <- Xwith %*% beta.old
    probab   <- offset * exp(eta)
    mu       <- probab
    V        <- mu
    deriv.mu <- offset * exp(eta)
    r.stand  <- (y - mu) / sqrt(V)

    if (weights.x){
      w.x <- sqrt(1 - hat(X))
    } else {
      w.x <- rep(1, length = nb)
    }

    g.objective <- function(beta) {
      eta      <- Xwith %*% beta
      probab   <- offset * exp(eta)
      mu       <- probab
      V        <- probab
      r.stand  <- (y - mu) / sqrt(V)
      deriv.mu <- offset * exp(eta)
      jinf     <- floor(mu - k * sqrt(V))
      jsup     <- floor(mu + k * sqrt(V))

      if (k == Inf){
        esp.cond <- rep(1, nb)
      }
      if (k != Inf){
        esp.cond <- - k * ppois(jinf, mu) + k * (1 - ppois(jsup, mu)) +
                    mu / sqrt(V) * (ppois(jinf, mu) - ppois(jinf - 1, mu) -
                    (ppois(jsup, mu) - ppois(jsup - 1, mu)))
      }

      a.const <- apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * esp.cond * deriv.mu),
                  2, sum)

      apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * basepsi(r.stand) * deriv.mu),
       2, sum) - a.const
    }

    grad.g <- function(beta) {
      delta <- .Machine$double.eps ^ .5
      Ident <- diag(1, length(beta))
      1 / delta * (apply(beta + delta * Ident, 2, g.objective) - as.vector(g.objective(beta)))
    }

    for (iiter in 1:100){
      g.old      <- g.objective(beta.old)
      grad.g.old <- grad.g(beta.old)
      csi        <- solve(grad.g.old, -g.old)
      beta.new   <- as.vector(beta.old + csi)
      convi      <- abs(max(beta.old - beta.new)) / abs(max(beta.old))

      done       <- (convi < mytol)
      if (done) { break }
      beta.old   <- beta.new
    }
    if (!done) {
      warning(paste("function for starting values used in mqpoisson failed to
                    converge in", 100, "steps."))
    }


    eta <- Xwith %*% beta.old
    fit <- offset * exp(eta)
    list(coef = beta.old, fitted.values = fit)
  }

  #Stopping rule
  irls.delta <- function(old, new) {
    abs(max(old - new)) / abs(max(old))
  }

  n <- length(case.weights)

  if (weights.x){
    w.x <- sqrt(1 - hat(x))
  } else {
    w.x <- rep(1, length = n)
  }

  #We fit the glm.rob for computing the starting values

  temp.rob <- glm.rob.poisson (X = x,
                               y = y,
                               weights.x = weights.x,
                               k = k,
                               offset = offset)

  resid.init <- y - temp.rob$fitted.values
  fit.init   <- temp.rob$fitted.values
  phi.init   <- 1

  done      <- FALSE
  conv      <- NULL
  qest      <- matrix(0, nrow = ncol(x), ncol = length(q))
  qfit      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar      <- matrix(0, nrow = ncol(x), ncol = length(q))
  huberised <- matrix(1, nrow = nrow(x), ncol = length(q))
  qphi      <- NULL

  for (i in 1:length(q)) {

    #We define the starting values
    resid <- resid.init
    fit   <- fit.init
    phi   <- phi.init
    w     <- case.weights
    coef  <- temp.rob$coef

    for (iiter in 1:maxit) {

      resid.old <- resid
      coef.old  <- coef

      # We define the expectation mu=t*exp(xb)
      probab   <- fit
      mu       <- probab
      deriv.mu <- mu

      #We define the variance
      V <- (phi ^ 2) * probab

      #We define the scale
      scale <- c(sqrt(V))

      #We standardize the residuals
      r.stand <- (y - mu) / sqrt(V)

      #we compute i1 and i2
      jinf <- floor(mu - k * sqrt(V))
      jsup <- floor(mu + k * sqrt(V))


      #We compute the values of a_j(b)
      if (k == Inf){
        a.j <- rep(1, n)
      }
      if (k != Inf){
        a.j <- (-k) * ppois(jinf, mu) + k * (1 - ppois(jsup, mu)) + mu / sqrt(V) *
               (ppois(jinf, mu) - ppois(jinf - 1, mu) - (ppois(jsup, mu) -
               ppois(jsup - 1, mu)))
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

      #we compute psi_q(r ) - E(psi_q(r ))
      A <- (hub.res - a.j)

      if (k == Inf){
        esp.carre.cond <- rep(1, n)
      }
      if (k != Inf){
        esp.carre.cond <- k * (ppois(jinf, mu) - ppois(jinf - 1, mu) +
                           (ppois(jsup, mu) - ppois(jsup - 1, mu))) +
                           (mu ^ 2 / V ^ (3 / 2)) * (ppois(jinf - 1, mu) -
                           ppois(jinf - 2 , mu) - (ppois(jinf, mu) -
                           ppois(jinf - 1, mu)) - (ppois(jsup - 1, mu) -
                           ppois(jsup - 2, mu)) + (ppois(jsup, mu) -
                           ppois(jsup - 1, mu))) + (mu / V ^ (3 / 2)) *
                           (ppois(jsup - 1, mu) - ppois(jinf, mu))
      }

      b.j <- 2 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
      B   <- diag(c(V * b.j))

      #We estimate betas
      temp <- coef + solve(t(x) %*% w %*% B %*% x) %*% t(x) %*% w %*% A
      coef <- temp

      eta   <- x %*% coef
      fit   <- offset * exp(eta)
      resid <- y - fit

      convi <- irls.delta(coef.old, coef)
      conv <- c(conv, convi)

      done <- (convi <= acc)
      if (done){
        break
      }

    }
    if (!done){
      warning(paste("mqpoisson failed to converge in", maxit, "steps at q = ", q[i]))
    }


    # Asymptotic estimated variance of the robust estimator
    probab   <- fit
    mu       <- probab
    deriv.mu <- mu

    #We define the variance
    V       <- phi * probab
    r.stand <- (y - mu) / sqrt(V)
    scale   <- c(sqrt(V))

    jinf <- floor(mu - k * sqrt(V))
    jsup <- floor(mu + k * sqrt(V))

    if (k == Inf) {
      esp.cond <- rep(1, n)
    } else {
      esp.cond <- -k * ppois(jinf, mu) + k * (1 - ppois(jsup, mu)) + mu / sqrt(V) *
                  (ppois(jinf, mu) - ppois(jinf - 1, mu) - (ppois(jsup, mu) -
                  ppois(jsup - 1, mu)))
    }

    esp.cond <- 2 * esp.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
    a.const  <- apply(x * as.vector(1 / n / sqrt(V) * w.x * esp.cond * deriv.mu), 2, sum)


    if (k == Inf) {
      esp.carre.cond <- rep(1, n)
    } else {
      esp.carre.cond <- k ^ 2 * (ppois(jinf, mu) + 1 - ppois(jsup, mu)) +
                         1 / V * (mu ^ 2 * (2 * ppois(jinf - 1, mu) -
                         ppois(jinf - 2, mu) - ppois(jinf, mu) - 2 *
                         ppois(jsup - 1, mu) + ppois(jsup - 2, mu) +
                         ppois(jsup, mu)) + mu * (ppois(jsup - 1, mu) -
                         ppois(jinf - 1, mu)))
    }

    esp.carre.cond <- 4 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) *
                       (r.stand <= 0)) ^ 2
    matQaux <- as.vector(esp.carre.cond / V * w.x ^ 2 * deriv.mu ^ 2)
    matQ1   <- (1 / n) * t(x) %*% (matQaux * x)
    matQ2   <- a.const %*% t(a.const)
    matQ    <- matQ1 - matQ2

    if (k == Inf) {
      esp.psi.score = 1 / sqrt(V)
    } else {
      esp.psi.score <- k * (ppois(jinf, mu) - ppois(jinf - 1, mu) +
                        (ppois(jsup, mu) - ppois(jsup - 1, mu))) +
                        (mu ^ 2 / V ^ (3 / 2)) * (ppois(jinf - 1, mu) -
                        ppois(jinf - 2, mu) - (ppois(jinf, mu) -
                        ppois(jinf - 1, mu)) - (ppois(jsup - 1, mu) -
                        ppois(jsup - 2, mu)) + (ppois(jsup, mu) -
                        ppois(jsup - 1, mu))) + (mu / V ^ (3 / 2)) *
                        (ppois(jsup - 1, mu) - ppois(jinf, mu))
    }

    esp.psi.score <- 2 * esp.psi.score  * (q[i] * (r.stand > 0) + (1 - q[i]) *
                      (r.stand <= 0))
    matMaux <- as.vector(esp.psi.score / sqrt(V) * w.x * deriv.mu ^ 2)
    matM    <- 1 / n * t(x) %*% (matMaux * x)
    matMinv <- solve(matM)

    as.var  <- 1 / n * matMinv %*% matQ %*% matMinv

    qest[, i] <- coef
    qfit[, i] <- fit
    qres[, i] <- y - fit
    qvar[, i] <- as.numeric(round(diag(as.var), 4))
  }

  if (compute.qscores) {
    q.scores <- mqpoisson.qscores(x = x,
                                  y = y,
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

  list(
    fitted.values = qfit,
    var.beta      = qvar,
    residuals     = qres,
    q.values      = q,
    coefficients  = qest,
    matQ          = matQ,
    matM          = matM,
    scale         = scale,
    q.scores      = q.scores,
    huberised.res = !(huberised)
    )
}

#' M-Quantile Estimation for Binary Data
#'
#' The function implements m-quantile estimation for binary data as proposed by Chambers et al. (2012).
#' The algorithm is based on robust quasi-likelihood estimation by Cantoni & Ronchetti (2001).
#' The function is called by the mquantreg-function,
#' when the method "binom" is selected. Important: Standalone use is possible, but not advised.
#'
#'@param x matrix, containing independent variables.
#'@param y vector, containing the dependent variable, usually a factor.
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'@param k a number greater than 0. k is the parameter for the huber psi-function, the default value is 1.6 for standalone use.
#'@param maxit an integer, which defines the maximal number of iteration
#'before the algorithm stops allthough it has not converged yet.
#'If the maximal number of iteration is reached a warning message will be prompted.
#'@param acc defines convergence criteria.
#'@param case.weights vector of observation weights;
#'if supplied, the algorithm fits to minimize the sum of the weights multiplied into the absolute residuals.
#'The length of weights must be the same as the number of observations.
#'The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous.
#'@param weights.x vector of weights, for estimation of robust glm to obtain initial weights in the discrete algorithms.
#'@param var.weights dataframe, if supplied the residuals are scaled.
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@param compute.qscores boolean, controls if q-scores are estimated.
#'@return See \code{\link{mquantreg}} for details.
#'\code{summary()}, \code{print()}, \code{fitted()} and \code{predict()}-methods are avaiable.
#'@references
#'\itemize{
#'\item Cantoni, E., & Ronchetti, E. (2001). Robust inference for generalized linear models. Journal of the American Statistical Association, 96(455), 1022-1030.
#'\item Chambers, R., Salvati, N., & Tzavidis, N. (2012). M-quantile regression for binary data with application to small area estimation.
#'}
#'
#'
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "binary",
#'                    measurement.error = 0.05)
#'
#'fit <- mquantreg(formula = Y ~ x1 + x2,data = df, q  = 0.5, method = "binom")
#'print(fit)
mqbinom <- function(x,
                    y,
                    q = 0.5,
                    k = 1.6,
                    maxit = 1000,
                    acc = 1e-04,
                    case.weights = NULL, #rep(1, nrow(x))
                    weights.x = FALSE,
                    var.weights = NULL, #rep(1, nrow(x))
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

    if(is.factor(y)) y <- as.numeric(y) - 1

    #We fit the glm.rob for computing the starting values
    temp.rob <- glm.rob.binom(x = x,
                              y = y,
                              weights.on.x = weights.x,
                              chuber = k
                              )

    resid.init <-  y - temp.rob$fitted.values
    fit.init   <- temp.rob$fitted.values

    done      <- FALSE
    conv      <- NULL
    qest      <- matrix(0, nrow = ncol(x), ncol = length(q))
    qfit      <- matrix(0, nrow = nrow(x), ncol = length(q))
    qres      <- matrix(0, nrow = nrow(x), ncol = length(q))
    qvar      <- matrix(0, nrow = ncol(x), ncol = length(q))
    huberised <- matrix(1, nrow = nrow(x), ncol = length(q))

    for(i in 1:length(q)) {

      #We define the starting values
      resid <- resid.init
      fit   <- fit.init
      # a.j   <- case.weights
      w     <- case.weights
      coef  <- temp.rob$coef

      for (iiter in 1:maxit) {
        resid.old <- resid
        coef.old  <- coef

        # We define the probability mu=exp(xb)/(1+exp(xb))
        probab <- fit
        mu     <- probab

        #We define the variance
        V <- probab * (1 - probab)

        #We define the scale
        scale <- c(sqrt(V))

        #We standardize the residuals
        r.stand <- (y - mu) / sqrt(V)

        #we compute i1 and i2
        jinf <- floor(mu - k * sqrt(V))
        jsup <- floor(mu + k * sqrt(V))
        ni   <- rep(1, n)

        #We compute the values of a_j(b)
        if (k == Inf) {
          a.j <- ni
        }
        if (k != Inf) {
          indic <- ifelse(jinf + 1 <= 1 & jsup >= 1, 1, 0)
          a.j   <- -k * pbinom(jinf, 1, probab) + k *
                    (1 - pbinom(pmin(jsup, ni), 1, probab)) + 1 / sqrt(V) *
                    ifelse(ni == 1, probab * indic, mu *
                    (pbinom(pmin(jsup - 1, ni - 1), pmax(0, 1), probab) -
                    pbinom(jinf - 1, pmax(0, 1), probab))) -
                    mu / sqrt(V) * (pbinom(pmin(jsup, 1), 1, probab) -
                    pbinom(jinf, 1, probab))
        }
        a.j <- 2 * a.j * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))

        #we define a part of w_j
        w <- c((mu * (1 - mu)) / scale) * c(w.x)

        #we compute psi_q(res)
        hub.res.sym        <- psi.huber((resid) / scale, k = k) * case.weights *
                               ((resid) / scale)
        huberised[, i]     <- (psi.huber((resid) / scale, k = k) == 1)
        hub.res            <- 2 * (1 - q[i]) * hub.res.sym
        hub.res[resid > 0] <- 2 * q[i] * hub.res.sym[resid > 0]


        #we compute psi_q(r )-E(psi_q(r ))
        A <- (hub.res - a.j)

        #We compute the -E(psi_prime)
        hub.res.sym        <- psi.huber((1 - mu) / scale, k = k) * case.weights *
                               ((1 - mu) / scale)
        hub.res            <- 2 * (1 - q[i]) * hub.res.sym
        hub.res[resid > 0] <- 2 * q[i] * hub.res.sym[resid > 0]
        B.hub.res          <- hub.res

        hub.res.sym        <- psi.huber(mu / scale, k = k) * case.weights *
                               (mu / scale)
        hub.res            <- 2 * (1 - q[i]) * hub.res.sym
        hub.res[resid > 0] <- 2 * q[i] * hub.res.sym[resid > 0]
        B                  <- c(c(V) * (B.hub.res + hub.res))

        #We estimate betas

        coef  <- coef + solve(t(x * w * B) %*% x) %*% t(x * w) %*% A
        eta   <- x %*% coef
        fit   <- exp(eta) / (1 + exp(eta))
        resid <- y - fit

        convi <- irls.delta(coef.old, coef)
        conv  <- c(conv, convi)

        done  <- (convi <= acc)

        if (done) { break }
      }
      if (!done) {
        warning(paste("mqbinom failed to converge in", maxit, "steps at q = ", q[i]))
      }



      # Asymptotic estimated variance of the robust estimator

      probab   <- fit
      mu       <- probab
      deriv.mu <- exp(eta) / ((1 + exp(eta)) ^ 2)

      #We define the variance
      V       <- probab * (1 - probab)
      r.stand <- (y - mu) / sqrt(V)

      basepsiq <- function(res, q) {
        hub.res.sym      <- psi.huber(res, k = k) * res
        hub.res          <- 2 * (1 - q) * hub.res.sym
        hub.res[res > 0] <- 2 * q * hub.res.sym[res > 0]
        return(hub.res)
      }

      # Asymptotic estimated variance of the robust estimator
      esp.cond       <- basepsiq((1 - mu) / sqrt(V), q = q[i]) * mu +
                          basepsiq(-mu / sqrt(V), q = q[i]) * (1 - mu)
      a.const        <- apply(x * as.vector(1 / n / sqrt(V) * esp.cond * deriv.mu), 2, sum)
      esp.carre.cond <- (basepsiq((1 - mu) / sqrt(V), q = q[i])) ^ 2 * mu +
                          (basepsiq(-mu / sqrt(V), q = q[i])) ^ 2 * (1 - mu)
      matQaux        <- as.vector(esp.carre.cond / V * w.x ^ 2 * (mu * (1 - mu))^2)
      matQ1          <- (1 / n) * t(x) %*% (matQaux * x)
      matQ2          <- a.const %*% t(a.const)
      matQ           <- matQ1 - matQ2

      esp.psi.score <- basepsiq((1 - mu) / sqrt(V), q = q[i]) - basepsiq(-mu / sqrt(V), q = q[i])
      matMaux       <- as.vector(esp.psi.score / sqrt(V) * w.x * V * (mu * (1 - mu)))
      matM          <- 1 / n * t(x) %*% (matMaux * x)
      matMinv       <- solve(matM)

      as.var        <- 1 / n * matMinv %*% matQ %*% matMinv


      qest[, i] <- coef
      qfit[, i] <- fit
      qres[, i] <- y - fit
      qvar[, i] <- as.numeric(round(diag(as.var), 4))
    }

    # compute q-values
    if (compute.qscores) {
      q.scores <- mqbinom.qscores(x,
                                  y,
                                  case.weights = case.weights,
                                  var.weights = var.weights,
                                  maxit = maxit,
                                  acc = acc,
                                  weights.x = weights.x,
                                  qgrid = qgrid,
                                  k.value = k)$qscores
    } else{
      q.scores <- NULL
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






# definition function glm.rob.binom
glm.rob.binom = function(x,
                         y,
                         weights.on.x = FALSE,
                         chuber = 1.6,
                         maxit = maxit) {

  y     <- as.numeric(as.character(y))
  mytol <- .Machine$double.eps ^ .25
  nb    <- length(y)

  basepsi <- function(x) {
    x * pmin(1, chuber / abs(x))
  }

  basepsiprime <- function(x) {
    1 * (abs(x) < chuber)
  }

  # Initialise
  beta.old <- as.vector(glm(y ~ x - 1, family = binomial)$coeff)
  Xwith    <- x
  eta      <- Xwith %*% beta.old
  probab   <- exp(eta) / (1 + exp(eta))
  mu       <- probab
  V        <- probab * (1 - probab)
  deriv.mu <- exp(eta) / ((1 + exp(eta)) ^ 2)
  r.stand  <- (y - mu) / sqrt(V)

  if (weights.on.x) {
    w.x <- sqrt(1 - hat(x))
  } else {
    w.x <- rep(1, length = nb)
  }

  g.objective <- function(beta) {
    eta     <- Xwith %*% beta
    probab  <- exp(eta) / (1 + exp(eta))
    mu      <- probab
    V       <- probab * (1 - probab)
    r.stand <- ifelse(sqrt(V) == 0, 0, (y - mu) / sqrt(V))

    deriv.mu <- exp(eta) / ((1 + exp(eta)) ^ 2)

    jinf <- floor(mu - chuber * sqrt(V))
    jsup <- floor(mu + chuber * sqrt(V))
    ni   <- rep(1, nb)

    if (chuber == Inf) {
      esp.cond <- numeric(nb)
    }
    if (chuber != Inf) {
      indic    <- ifelse(jinf + 1 <= 1 & jsup >= 1, 1, 0)
      esp.cond <- -chuber * pbinom(jinf, 1, probab) + chuber *
        (1 - pbinom(pmin(jsup, ni), 1, probab)) +
        1 / sqrt(V) * ifelse(ni == 1, probab * indic, mu *
                               (pbinom(pmin(jsup - 1, ni - 1), pmax(0, 1), probab) -
                                  pbinom(jinf - 1, pmax(0, 1), probab))) - mu / sqrt(V) *
        (pbinom(pmin(jsup, 1), 1, probab) - pbinom(jinf, 1, probab))
    }

    a.const <- apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * esp.cond * deriv.mu),
                     2, sum)

    apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * basepsi(r.stand) * deriv.mu),
          2 ,sum) - a.const
  }

  grad.g <- function(beta) {
    delta <- .Machine$double.eps ^ .5
    Ident <- diag(1, length(beta))
    1 / delta * (apply(beta + delta * Ident, 2, g.objective) -
                   as.vector(g.objective(beta)))
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
    warning(paste("function for starting values used in mqbinom failed to converge
                  in", 100, "steps."))
  }

  eta <- Xwith %*% beta.new
  fit <- exp(eta) / (1 + exp(eta))
  list(coef = beta.new, fitted.values = fit)
}


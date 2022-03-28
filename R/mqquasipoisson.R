#' M-Quantile Estimation for Count Data with Quasi Poisson
#'
#' The algorithm is based on robust quasi-likelihood estimation by Cantoni & Ronchetti (2001).
#' The function is called by the mquantreg-function, when the method "quasipoisson" is selected.
#' Important: Standalone use is possible, but not advised.
#'
#'@param x matrix, containing independent variables.
#'@param y vector, containing dependent variable.
#'@param q number strictly between 0 and 1, which specifies the m-quantile to be estimated.
#'@param k a number greater than 0. k is the parameter for the huber psi-function, the default value is 1.345.
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
#'@param phi.init numeric
#'@param qgrid vector of quantiles, which are to be used for qscore calculations.
#'@param compute.qscores boolean, controls if q-scores are estimated.
#'@return See \code{\link{mquantreg}} for details.
#'@references
#'Cantoni, E., & Ronchetti, E. (2001). Robust inference for generalized linear models. Journal of the American Statistical Association, 96(455), 1022-1030.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "count.poisson",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2 * x1",data = df, q  = 0.5, method = "quasipoisson")
#'print(fit)
mqquasipoisson <- function(x,
                           y,
                           q = 0.5,
                           k = 1.6,
                           maxit,
                           acc = 1e-04,
                           offset = NULL,
                           case.weights = NULL,
                           weights.x = FALSE,
                           var.weights = NULL,
                           phi.init = 1,
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
                           ){

  if (is.null(offset)) {
    offset <- rep(1, nrow(x))
  }

  if (is.null(case.weights)){
    case.weights <- rep(1, nrow(x))
  }


  if (is.null(var.weights)) {
    var.weights <- rep(1, nrow(x))
  }

  Ez <- function(q, k){
    N      <- 10000000
    e      <-  rnorm(N, 0, 1)
    res.q  <- 2 * psi.huber(e, k = k) * e * (q * (e > 0) + (1 - q) * (e <= 0))
    mean(res.q ^ 2)
  }

  glm.rob.poisson <- function(x,
                              y,
                              weights.on.x = FALSE,
                              chuber = 1.345,
                              offset = offset){

      mytol  <-  .Machine$double.eps ^ .25
      nb     <-  length(y)


      basepsi <- function(x){
        x * pmin(1, chuber / abs(x))
      }

      basepsiprime <- function(x){
        1 * (abs(x) < chuber)
      }

      # Initializing....
      beta.old <- as.vector(glm(y ~ x - 1, family = poisson, offset = log(offset))$coeff)
      Xwith    <- x
      eta      <- Xwith %*% beta.old
      probab   <- offset * exp(eta)
      mu       <- probab
      V        <- mu
      deriv.mu <- offset * exp(eta)
      r.stand  <- (y - mu) / sqrt(V)

      if (weights.on.x) {
        w.x <- sqrt(1 - hat(x))
      } else {
        w.x <- rep(1, length = nb)
      }

      g.objective <- function(beta){
        eta      <- Xwith %*% beta
        probab   <- offset * exp(eta)
        mu       <- probab
        V        <- probab
        r.stand  <- (y - mu) / sqrt(V)
        deriv.mu <- offset * exp(eta)
        jinf     <- floor(mu - chuber * sqrt(V))
        jsup     <- floor(mu + chuber * sqrt(V))

        if (chuber == Inf) {
          esp.cond <- rep(1, nb)
        }
        if (chuber != Inf) {
          esp.cond <- -chuber * ppois(jinf, mu) + chuber * (1 - ppois(jsup, mu)) +
                      mu / sqrt(V) * (ppois(jinf, mu) - ppois(jinf - 1, mu) -
                      (ppois(jsup, mu) - ppois(jsup - 1, mu)))
        }
        a.const <- apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * esp.cond * deriv.mu), 2, sum)

        apply(Xwith * as.vector(1 / nb / sqrt(V) * w.x * basepsi(r.stand) * deriv.mu), 2, sum) - a.const
      }

      grad.g = function(beta){
        delta <- .Machine$double.eps ^ .5
        Ident <- diag(1, length(beta))

        1 / delta * (apply(beta + delta * Ident, 2, g.objective) - as.vector(g.objective(beta)))
      }

      iter.cnt <- 0
      # Main
      repeat {
        iter.cnt   <- iter.cnt + 1
        g.old      <- g.objective(beta.old)
        grad.g.old <- grad.g(beta.old)
        csi        <- solve(grad.g.old, -g.old)
        beta.new   <- as.vector(beta.old + csi)

        if (abs(max(beta.old - beta.new)) / abs(max(beta.old)) < mytol)
          break
        if (iter.cnt > 100)
          break
        beta.old = beta.new
        NULL
      }

      eta <- Xwith %*% beta.old
      fit <- offset * exp(eta)

      list(coef = beta.old, fitted.values = fit)
    }

  #Stopping rule
  irls.delta <- function(old, new){
    abs(max(old - new)) / abs(max(old))
  }


  n <- length(case.weights)
  # print(n)
  if (weights.x) {
    w.x <- sqrt(1 - hat(x))
  } else {
    w.x <- rep(1, length = n)
  }

  #We fit the glm.rob for computing the starting values
  p <- ncol(x)

  temp.rob <- glm.rob.poisson(x = x,
                              y = y,
                              weights.on.x = weights.x,
                              chuber = k,
                              offset = offset)

  resid.init <- y - temp.rob$fitted.values
  fit.init   <- temp.rob$fitted.values

  done      <- FALSE
  conv      <- NULL
  qest      <- matrix(0, nrow = ncol(x), ncol = length(q))
  qfit      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres      <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar      <- matrix(0, nrow = ncol(x), ncol = length(q))
  qphi      <- NULL
  huberised <- matrix(1, nrow = nrow(x), ncol = length(q))

  for (i in 1:length(q)) {

    #We define the starting values
    resid <- resid.init
    fit   <- fit.init
    phi   <- phi.init
    a.j   <- case.weights
    w     <- case.weights
    coef  <- temp.rob$coef

    if (q[i] == 0.5) {
      ExpectedZ <- 2 * pnorm(k) - 1 - 2 * k * dnorm(k) + 2 * (k ^ 2) * (1 - pnorm(k))
    }

    if (q[i] != 0.5) {
      ExpectedZ <- Ez(q[i], k = k)
    }

    for (iiter1 in 1:maxit) {
      for (iiter in 1:maxit) {
        resid.old <- resid
        coef.old  <- coef


        # We define the probability mu=t*exp(xb)
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
        if (k == Inf) {
          a.j <- rep(1, n)
        }

        if (k != Inf) {
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

        #we compute psi_q(r )-E(psi_q(r ))
        A <- (hub.res - a.j)

        if (k == Inf){
          esp.carre.cond <- rep(1, n)
        }

        if (k != Inf){
          esp.carre.cond <- k * (ppois(jinf, mu) - ppois(jinf - 1, mu) +
                            (ppois(jsup, mu) - ppois(jsup - 1, mu))) + (mu ^ 2 / V ^ (3 / 2)) *
                            (ppois(jinf - 1, mu) - ppois(jinf - 2, mu) - (ppois(jinf, mu) -
                            ppois(jinf - 1, mu)) - (ppois(jsup - 1, mu) - ppois(jsup - 2, mu)) +
                            (ppois(jsup, mu) - ppois(jsup - 1, mu))) + (mu / V ^ (3 / 2)) *
                            (ppois(jsup - 1, mu) - ppois(jinf, mu))
        }
        b.j <- 2 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
        B   <- diag(c(V * b.j))

        #We estimate betas
        temp  <- coef + solve(t(x) %*% w %*% B %*% x) %*% t(x) %*% w %*% A
        coef  <- temp
        eta   <- x %*% coef
        fit   <- offset * exp(eta)
        resid <- y - fit

        convi <- irls.delta(coef.old, coef)
        conv  <- c(conv, convi)
        done  <- (convi <= acc)
        if (done)
          break
      }

      #estimation of phi

      basepsiq <- function(res, q){
        hub.res.sym      <- psi.huber(res, k = k) * res
        hub.res          <- 2 * (1 - q) * hub.res.sym
        hub.res[res > 0] <- 2 * q * hub.res.sym[res > 0]
        return(hub.res)
      }

      diff    <- 1
      phi.old <- phi

      while (diff > acc){
        res.std <- (y - fit) / (phi.old * sqrt(fit))
        res.q   <- basepsiq(res = res.std, q = q[i])
        wi0     <- res.q / res.std

        phi.new <- (1 / (n - p)) * sum((wi0 ^ 2 * ((y - fit) / sqrt(fit)) ^ 2) / ExpectedZ)
        diff    <- (phi.new - phi.old) ^ 2
        phi.old <- phi.new
      }

      convi1 <- (phi - phi.new) ^ 2
      phi    <- phi.new
      done1  <- (convi1 <= acc)
      if (done1)
        break
    }
    if (!done)
      warning(paste("MQQuasi-Poisson beta failed to converge in", maxit, "steps at q = ", q[i]))
    if (!done1)
      warning(paste("MQQuasi- Poisson phi failed to converge in", maxit, "steps at q = ", q[i]))

    # # Asymptotic estimated variance of the robust estimator
    probab   <- fit
    mu       <- probab
    deriv.mu <- mu

    #We define the variance
    V       <- (phi ^ 2) * probab
    r.stand <- (y - mu) / sqrt(V)
    scale   <- c(sqrt(V))

    jinf <- floor(mu - k * sqrt(V))
    jsup <- floor(mu + k * sqrt(V))

    if (k == Inf){
      esp.cond <- rep(1, n)
    } else {
      esp.cond <- -k * ppois(jinf, mu) + k * (1 - ppois(jsup, mu)) + mu / sqrt(V) *
        (ppois(jinf, mu) - ppois(jinf - 1, mu) - (ppois(jsup, mu) -
                                                    ppois(jsup - 1, mu)))
    }

    esp.cond <- 2 * esp.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
    a.const  <- apply(x * as.vector(1 / n / sqrt(V) * w.x * esp.cond * deriv.mu), 2, sum)


    if (k == Inf) {
      esp.carre.cond <- 1
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

    if (k == Inf){
      esp.psi.score <- 1 / sqrt(V)
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
    qphi[i]   <- phi
    qvar[, i] <- as.numeric(round(diag(as.var), 4))
  }



  if (compute.qscores) {
    q.scores <- mqquasipoisson.qscores(x,
                                       y,
                                       offset = offset,
                                       case.weights = case.weights,
                                       maxit = maxit,
                                       acc = acc,
                                       weights.x = weights.x,
                                       qgrid = qgrid,
                                       k.value = k,
                                       epsilon = 0.01)$qscores
  } else {
    q.scores <- NULL
  }

  list(
    fitted.values = qfit,
    var.beta      = qvar,
    residuals     = qres,
    q.values      = q,
    coefficients  = qest,
    q.phi         = qphi,
    scale         = scale,
    q.scores      = q.scores,
    huberised.res = !(huberised)
  )
}

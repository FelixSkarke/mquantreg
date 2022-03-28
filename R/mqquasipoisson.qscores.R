mqquasipoisson.qscores <- function(x,
                                   y,
                                   offset = rep(1, nrow(x)),
                                   case.weights = rep(1, nrow(x)),
                                   maxit = 1000,
                                   acc = 1e-4,
                                   weights.x = FALSE,
                                   qgrid = c(0.1,
                                             0.15,
                                             0.2,
                                             0.25,
                                             0.3,
                                             0.35,
                                             0.4,
                                             seq(from = 0.45, to = 0.55, by = 0.005),
                                             0.6,
                                             0.65,
                                             0.7,
                                             0.75,
                                             0.8,
                                             0.85,
                                             0.9
                                             ),
                                   k.value = 1.6,
                                   epsilon = 0.01){

    temp <- mqquasipoisson(x = x,
                           y = y,
                           offset = offset,
                           case.weights = case.weights,
                           maxit = maxit,
                           acc = acc,
                           q = qgrid,
                           k = k.value,
                           weights.x = weights.x)

    qfit <- temp$fitted.values
    qres <- temp$residuals
    qest <- temp$coefficients
    lfit <- x %*% qest
    Q    <- length(qgrid)

    midgrid      <- round((Q + 1) / 2)
    lhalf        <- lfit[, midgrid]
    lhalf.order  <- order(lhalf)
    qhalf        <- qfit[, midgrid]
    ystar        <- NULL
    zero         <- sum(y == 0)
    ystar[y > 0] <- log(y[y > 0]) - log(offset[y > 0])

    if (zero > 0){
      tmp <- which(y == 0)
      for (tmp1 in tmp) {
        ystar[tmp1] <- min((log(1 - epsilon) + log(offset[tmp1])),-lhalf[tmp1]) -
                       2 * log(offset[tmp1])
      }
    }


    lqvals     <- main.qscores(y = ystar, yhatq = lfit, qvals = qgrid)
    lqvals.std <- ((lqvals * 0.5) / mean(lqvals))

    list(
      fitted.values = qfit,
      residuals     = qres,
      coefficients  = qest,
      qscores       = lqvals,
      qscores.std   = lqvals.std
      )
  }

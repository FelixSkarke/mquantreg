mqbinom.qscores <- function(x,
                            y,
                            case.weights = rep(1, nrow(x)),
                            var.weights = rep(1, nrow(x)),
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
                            plots = FALSE,
                            k.value = 1.6) {

    temp <- mqbinom(y = y,
                    x = x,
                    case.weights = case.weights,
                    var.weights = var.weights,
                    maxit = maxit,
                    acc = acc,
                    q = qgrid,
                    k = k.value,
                    weights.x = weights.x
                    )

    qfit <- temp$fitted.values
    qres <- temp$residuals
    qest <- temp$coefficients
    lfit <- x %*% qest
    Q    <- length(qgrid)

    midgrid     <- round((Q + 1) / 2)
    lhalf       <- lfit[, midgrid]
    lhalf.order <- order(lhalf)
    qhalf       <- qfit[, midgrid]
    ystar       <- log((0.5 * (y + qhalf)) / (1 - 0.5 * (y + qhalf)))
    lqvals      <- main.qscores(y = ystar, yhatq = lfit, qvals = qgrid)
    lqvals.std  <- lqvals * 0.5 / mean(lqvals)

    if (plots) {
      #quartz(width=15,height=5)

      par(mfrow = c(1, 3))
      plot(x = range(lhalf),
           y = c(-0.1, 1.1),
           type = "n",
           main = "QLogistic Fits",
           xlab = "Linear Component of Logistic Fit",
           ylab = "Y"
           )
      points(x = lhalf[y == 1], y = 0.1 + y[y == 1], col = "blue")
      points(x = lhalf[y == 0],
             y = y[y == 0] - 0.1,
             col = "green")
      for (q in 1:Q)
        lines(x = lhalf[lhalf.order], y = qfit[lhalf.order, q])
        lines(x = lhalf[lhalf.order],
              y = qhalf[lhalf.order],
              col = "red",
              lwd = 2
              )
      abline(h = 0.5, lty = 3)
      boxplot(list("Y=1" = lqvals[y == 1], "Y=0" = lqvals[y == 0]), main = "QLogistic Coefficients")
      abline(h = 0.5, lty = 3)
      plot(x = range(lhalf),
           y = range(lqvals),
           type = "n",
           main = "QLogistic Coefficients",
           xlab = "Linear Component of Logistic Fit",
           ylab = ""
           )
      points(x = lhalf[y == 1], y = lqvals[y == 1], col = "blue")
      points(x = lhalf[y == 0], y = lqvals[y == 0], col = "green")
      abline(h = 0.5, lty = 3)
      par(mfrow = c(1, 1))
    }

    list(fitted.values = qfit,
         residuals = qres,
         coefficients = qest,
         qscores = lqvals,
         qscores.std = lqvals.std
         )
  }

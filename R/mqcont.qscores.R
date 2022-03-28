mqcont.qscores <- function(x,
                           y,
                           offset = rep(1, nrow(x)),
                           case.weights = rep(1, nrow(x)),
                           var.weights = rep(1, nrow(x)),
                           maxit = 1000,
                           acc = 1e-4,
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
                                     0.9
                                     ),
                            k.value = 1.345,
                            epsilon = 0.01,
                            scale.estimator = 'Mad'){
  temp = mqcont(x = x,
                y = y,
                offset = offset,
                case.weights = case.weights,
                maxit = maxit,
                acc = acc,
                q = qgrid,
                k = k.value,
                scale.estimator = scale.estimator
                 )

  qfit   <- temp$fitted.values
  qres   <- temp$residuals
  qest   <- temp$coefficients
  lfit   <- x %*% qest
  lqvals <- main.qscores(y = y,
                         yhatq = lfit,
                         qvals = qgrid)
  lqvals.std <- ((lqvals * 0.5) / mean(lqvals))

  return(list(fitted.values = qfit,
              residuals     = qres,
              coefficients  = qest,
              qscores       = lqvals,
              qscores.std   = lqvals.std
  )
  )
}

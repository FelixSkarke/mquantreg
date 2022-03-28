#' Fitted vs. Residuals Plot
#'
#' The fit.res plot creates a residuals over fitted value plot.
#' The smoothed residual indicates the mean residual conditional on the fitted values.
#' Additionally, observation which are downweighted are marked with a triangle.
#'
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param sqr.res boolean controlling if the plot should show the squared residuals instead.
#'@param add.qvals numeric, containing additional m-quantiles to be shown in plot.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data( n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                   measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2 * x1",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "fit.res")
#'
#'@export
ggplot.fit.res <- function(fit,
                           sqr.res = FALSE,
                           add.qvals = FALSE) {
  # Adding q-values if extra a passed and estimating models with the extra q-values
  qvals <- fit$q.values
  if (all(add.qvals != FALSE)) {
    qvals <- c(add.qvals, qvals)
    fit <- mquantreg(formula = fit$formula,
                     data = fit$data,
                     q  = qvals,
                     method = fit$method,
                     k = fit$k,
                     scale.estimator = fit$scale.estimator
                     )
    }

  # Sorting q-values
  qvals <- sort(qvals)

  # Creating dataframe for plotting
  data_df <- data.frame()
  for (i in 1:length(qvals)) {
    data_df <- rbind(data_df,
                     cbind(fit$fitted.values[, i],
                           fit$residuals[, i],
                           fit$huberised.res[, i],
                           rep(qvals[i], length(fit$residuals[, i]))
                     ))
  }

  colnames(data_df) <- c("fitted", "residuals", "huberised", "q")
  data_df$q         <- as.factor(data_df$q)
  data_df$huberised <- as.factor(as.logical(data_df$huberised))

  # Plot with squared residuals
  if (sqr.res) {
    data_df$residuals.sqr <- data_df$residuals ^ 2

    ggplot(data_df, aes(residuals.sqr, fitted)) +
      geom_point(aes(col = q, shape = huberised)) +
      geom_smooth(aes(col = q), se = FALSE, method = loess)
  } else {
    # Plot with residuals
    ggplot(data_df, aes(residuals, fitted, col = q)) +
      geom_point(aes( shape = huberised)) +
      geom_smooth(se = FALSE,
                  method = loess,
                  size = 1.2) +
      xlab("fitted values") +
      ylab("residuals") +
      labs(color = "q-value")
  }
}





#' Proportion of Residuals smaller than 0 over the Fitted Values
#'
#' This plot shows the proportion of residuals smaller than 0 over the fitted values.
#' This can give an indication of model fit and make m-quantile regression comparable with quantile regression.
#'
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param add.qvals numeric, containing additional m-quantiles to be shown in plot.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                   measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "prop.res")
#'
#'@export
ggplot.prop.res <- function(fit,
                            add.qvals = FALSE) {
  # Adding q-values if extra a passed and sorting q-values
  qvals <- fit$q.values
  if (all(add.qvals != FALSE)) {
    qvals <- c(add.qvals, qvals)

    # Sorting q-values
    qvals <- sort(qvals)
    fit   <- mquantreg(formula = fit$formula,
                     data = fit$data,
                     q  = qvals,
                     method = fit$method,
                     k = fit$k,
                     scale.estimator = fit$scale.estimator
                     )
  }



  # Creat dataframe for plotting
  resp   <- fit$residuals < 0
  data_df<- data.frame()

  for (i in 1:length(qvals)) {
    glm.fit    <- glm(resp[, i] ~ fit$X[, -1], family = "binomial")
    mq.fitted  <- fit$fitted.values[, i]
    glm.fitted <- glm.fit$fitted.values

    data_df <- rbind(data_df, cbind(mq.fitted, glm.fitted, rep(qvals[i], length(fit$residuals[, i]))))
  }

  colnames(data_df) <- c("mq.fitted", "glm.fitted", "q")
  data_df$q         <- as.factor(data_df$q)
  rownames(data_df) <- seq(length = nrow(data_df))

  # Plot
  ggplot(data_df, aes(mq.fitted, glm.fitted, col = q)) +
    geom_line(size = 1.2) +
    ylim(c(0, 1)) +
    xlab("fitted values") +
    ylab(expression(paste("Proportion of ", r[q], " < 0"))) +
    labs(color = "q-value")
}


#' Fitted vs. Observed Values
#'
#' This plot shows the fitted values over the observed. This gives a good indication of how good model fit is-
#' Additionally observation which are downweighted are marked with a triangle.
#'
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param add.qvals numeric, containing additional m-quantiles to be shown in plot.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data( n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2 * x1",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "fitted.observed")
#'@export
ggplot.fitted.observed <- function(fit,
                                   add.qvals = FALSE) {
  fit <- fit2
  add.qvals <- 0.4
  # Adding q-values if extra a passed and estimating models with the extra q-values
  qvals <- fit$q.values
  if (all(add.qvals != FALSE)) {
    qvals <- c(add.qvals, qvals)
    fit <- mquantreg(formula = fit$formula,
                     data = fit$data,
                     q  = qvals,
                     method = fit$method,
                     k = fit$k,
                     scale.estimator = fit$scale.estimator
                     )
  }

  # Sorting q-values
  qvals <- sort(qvals)

  # Creating dataframe for plotting
  data_df <- data.frame()
  for (i in 1:length(qvals)) {
    data_df <- rbind(data_df,
                     cbind(fit$fitted.values[, i],
                           fit$Y,
                           fit$huberised.res[, i],
                           rep(qvals[i], length(fit$residuals[, i]))
                           )
                     )
  }

  colnames(data_df) <- c("fitted", "true.y", "huberised", "q")
  data_df$q         <- as.factor(data_df$q)
  data_df$huberised <- as.factor(as.logical(data_df$huberised))

  # Plot
  ggplot(data_df, aes(fitted, true.y)) +
    geom_abline(
      slope = 1,
      size = 1.0,
      alpha = 0.5,
      linetype = 2
    ) +
    geom_point(aes(col = q, shape = huberised)) +
    geom_smooth(aes(col = q),
                se = FALSE,
                method = loess,
                size = 1.2) +
    xlab("fitted values") +
    ylab("obseverved values") +
    labs(color = "q-value") +
    ylim(c(min(data_df$fitted, data_df$true.y),
           max(data_df$fitted, data_df$true.y)
           )) +
    xlim(c(min(data_df$fitted, data_df$true.y),
           max(data_df$fitted, data_df$true.y)
           ))
}


#' Optimal k-Plot
#'
#' This plot shows the distance of the residual distribution to the normal distribution.
#' It there by uses the method proposed by James Dawber (2017). The optimal k is following
#' the method the k with the smallest distance, hence the minimum in the plot.
#'
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param kstart numeric, the first k estimated and shown in the plot.
#'@param kend  numeric, the last k estimated and shown in the plot.
#'@param kstep numeric, controls the steps in which k is increased between kstart and kend.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data(n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2 * x1",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "optk")
#'@export
ggplot.optk <- function(fit,
                        kstart = 0.5,
                        kend = 3,
                        kstep = 0.1) {
  ### Inverse M-quantile funciton for normal distribution
  Gms.cMAD <- function(y, s, k = 1.345) {
    numerator <- -k * pnorm(y - k * s) + (1 / s) * (-dnorm(y) + dnorm(y - k * s) -
                        y * (pnorm(y) - pnorm(y - k * s)))

    denominator <- numerator - ((1 / s) * (-dnorm(y + k * s) + dnorm(y) -
                          y * (pnorm(y + k * s) - pnorm(y))) + k - k *
                          pnorm(y + k * s))
    return(numerator / denominator)
  }

  kgrid <- seq(kstart, kend, kstep)
  qvals <- seq(0.01, 0.99, 0.01)
  k     <- fit$k
  fit   <- mquantreg(formula = fit$formula,
                     data = fit$data,
                     q  = qvals,
                     method = fit$method,
                     k = fit$k,
                     scale.estimator = fit$scale.estimator
                     )

  data_df <- data.frame()
  for (i in 1:length(kgrid)) {
    fit_k_i <- mquantreg(formula = fit$formula,
                         data = fit$data,
                         q  = qvals,
                         method = fit$method,
                         k = kgrid[i],
                         scale.estimator = fit$scale.estimator
                         )

    fit_k_50 <- mquantreg(formula = fit$formula,
                          data = fit$data,
                          q  = 0.5,
                          method = fit$method,
                          k = kgrid[i],
                          scale.estimator = fit$scale.estimator
                          )

    d_j <- sum((Gms.cMAD(apply(fit_k_i$residuals, 2, median) -
              median(fit_k_50$residuals), fit_k_i$scale, k = 3) - qvals) ^ 2)

    data_df <- rbind(data_df, cbind(kgrid[i], d_j))
  }

  colnames(data_df) <- c("k", "d_j")

  # Plot
  ggplot(data_df, aes(k, d_j)) +
    geom_line(aes(group = 1), size = 1.2) +
    xlab("k-values") +
    ylab(expression(d[j])) +
    geom_vline(xintercept = k, linetype = 2) +
    annotate(
      "text",
      x = k + 0.05,
      y = min(d_j) * 0.995,
      label = paste("k = ", k),
      hjust = 0
    )

}

#' Coefficients over the M-Quantiles
#'
#' This plot shows the coefficients over the m-quantiles.
#' This can show the changing relationship for changing q.
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param coef numeric, coefficient to plot. Defaults to \code{FALSE}.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data( n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                   measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2 * x1",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "coef.q", coef = 0)
#'@export
ggplot.coef.q <- function(fit , coef = FALSE) {
  fit <- fit
  coef <- FALSE
  qvals <- seq(0.01, 0.99, 0.01)

  if (coef == FALSE){
    coef <- 2:(ncol(fit$coefficients))
  } else { coef <- coef + 1 }

  fit_qval <- mquantreg(formula = fit$formula,
                        data = fit$data,
                        q  = qvals,
                        method = fit$method,
                        scale.estimator = fit$scale.estimator
                        )

  # Creating dataframe for plotting
  data_df <- data.frame()


  for (i in 1:length(qvals)) {
    for (coef.i in 1:length(coef)) {
      j <- coef[coef.i]
      # ,beta_exp = NA,
      data_ij <- data.frame(coef_div = NA,
                            j = NA,
                            q = NA
                            )
      data_ij$coef_div <- as.numeric(((fit_qval$coefficients[i, j]) - mean(fit_qval$coefficients[, j])) /
                                       mean(fit_qval$coefficients[, j]))  * 100
      data_ij$j        <- j
      # data_ij$beta_exp <- expression(paste0("beta[", j - 1, "]"))
      data_ij$q        <- as.numeric(qvals[i])
      data_df          <- rbind(data_df, data_ij)
    }
  }

  data_df$j        <- as.factor(data_df$j)
  data_df$coef_div <- as.numeric(as.character(data_df$coef_div))
  data_df$q        <- as.numeric(as.character(data_df$q))

  ggplot(data_df, aes(q, coef_div, color = j)) +
    geom_line(size = 1.2) +
    xlab("q") +
    ylab("% divergence from mean coefficient") +
    xlim(c(0, 1)) +
    labs(color = "coefficient")
}

#' Proportion of residuals huberised over k
#'
#' This plot shows the effect of k on the amount of observation being classified as outliers and therefore downweighted.
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@param kstart numeric, the first k estimated and shown in the plot.
#'@param kend numeric, the last k estimated and shown in the plot.
#'@param kstep numeric, controls the steps in which k is increased between kstart and kend.
#'@param add.qvals numeric, containing additional m-quantiles to be shown in plot.
#'
#'@examples
#'library(mq1)
#'
#'df <- simulate_data( n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "prop.huber.k")
#'@export
ggplot.prop.huber.k <- function(fit,
                                kstart = 0.5,
                                kend = 3,
                                kstep = 0.1,
                                add.qvals = FALSE) {
  kvals <- seq(kstart, kend, kstep)
  # Adding q-values if extra a passed and estimating models with the extra q-values
  qvals <- fit$q.values
  if (all(add.qvals != FALSE)) {
    qvals <- c(add.qvals, qvals)
  }

  # Creating dataframe for plotting
  data_df <- data.frame()
  for (i in 1:length(kvals)) {
    for (j in 1:length(qvals)) {
      fit_k <- mquantreg(formula = fit$formula,
                         data = fit$data,
                         q  = qvals[j],
                         method = fit$method,
                         k = kvals[i],
                         scale.estimator = fit$scale.estimator
                         )
      huber.props <- colMeans(fit_k$huberised.res)
      data_df <- rbind(data_df, cbind(huber.props, kvals[i], qvals[j]))
    }
  }

  colnames(data_df) <- c("huber.prop", "k", "q")
  data_df$q <- as.factor(data_df$q)

  # Plot
  ggplot(data_df, aes(k, huber.prop, color = q)) +
    geom_line(size = 1.2) +
    xlab("k") +
    ylab("Proportion of huberised residuals") +
    labs(color = "q-value") +
    ylim(c(0, 1)) +
    xlim(c(kstart, kend))
}

#' Proportion of residuals huberised over q
#'
#' This plot shows the effect of k on the amount of observation being classified as outliers and therefore downweighted.
#'@param fit a mquantfit object, produced by calling the mquantreg-function.
#'@examples
#'library(mq1)
#'
#'df <- simulate_data( n = 100,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'
#'fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
#'plot(fit, plottype = "prop.huber.q")
#'@export
ggplot.prop.huber.q = function(fit) {
  qvals <- seq(0.01, 0.99, 0.01)
  fit_q <- mquantreg(formula = fit$formula,
                     data = fit$data,
                     q  = qvals,
                     method = fit$method,
                     k = fit$k,
                     scale.estimator = fit$scale.estimator
                     )
  huber.props <- colMeans(fit_q$huberised.res)

  # Creating dataframe for plotting
  data_df <- data.frame()
  for (i in 1:length(qvals)) {
    data_df <- rbind(data_df, cbind(huber.props[i], qvals[i]))
  }

  colnames(data_df) <- c("huber.prop", "q")

  # Plot
  ggplot(data_df, aes(q, huber.prop)) +
    geom_line(size = 1.2) +
    xlab("q") +
    ylab("Proportion of huberised residuals") +
    labs(color = "q-value")  +
    ylim(c(0, 1)) +
    xlim(c(0, 1))
}

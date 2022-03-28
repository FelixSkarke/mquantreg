#' Wald test for M-quantile regression coefficients and linear hypotheses
#'
#' This function allows to test the Null hypothesis of regression coefficients being zero.
#' It only works for coefficients estimated with a likelihood approach for a continuous dependent variable.
#'
#' @param object an object of class "mquantfit", usually, a result of a call to \link{mquantreg}
#' @param id.beta identifier for the coefficient(s) to be jointly tested. Expects vector of positions
#' in coefficients of \code{object}. Defaults to "all", which does not include the intercept.
#' @param q identifies the m-quantile(s) for which coefficients should be tested.
#' Expects vector of m-quantiles e.g. \code{c(0.3, 0.5, 0.7)}.Defaults to "all".
#' @param R an optional matrix conformable to the coeffiecients from \code{object} in a way that
#' \code{R \%*\% beta.q} gives linear combinations of coefficients to be tested. Defaults to NULL.
#' @param H0 a vector specifying the null hypotheses.
#' @param each gives the user the possibility to test each coefficient individually. If set to \code{TRUE},
#' all chosen coefficients in id.beta will be tested.Defaults to FALSE.
#'
#' @return the function mqwald returns the the result of
#'
#' @param coefficents the estimated coeffients
#' @param varCov the variance-covariance matrix of the coefficients
#' @param statistic the test statistic
#' @param p.values the p-values for the Wald tests
#' @references
#' Bianchi, A., Fabrizi, E., Salvati, N., & Tzavidis, N. (2015).
#' M-quantile Regression: Diagnostics and the parametric representation of the model.
#' Book of Abstracts SIS - CLADAG, 303-306
#' Bianchi, A., & Salvati, N. (2015). Asymptotic properties and variance estimators
#' of the M-quantile regression coefficients estimators. Commun.Stat. Theory., 44, 2016–2429
#' @examples
#' df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                   measurement.error = 0.05)
#'
#' fit <- mquantreg(formula = Y ~ x1 + x2, data = df, q  = 0.5, method = "continuous")
#'
#' hyp_test <- mqwald(fit)
#' print(hyp_test)
#'
#' @export
mqwald <- function(object, id.beta = "all", q = "all", R = NULL, H0 = NULL, each = FALSE){

  cl <- match.call()
  print_id <- list("id.beta" = id.beta,
                   "q" = q,
                   "R" = R,
                   "H0" = H0,
                   "each" = each)

  if (object$method != "continuous") {
    stop("Wald hypothesis testing is only implemented for continuous dependent variables
         estimated with the likelihood method.")
  }

  if (!is.numeric(id.beta) & id.beta[1] != "all") {
    stop("id.beta must be an vector indicating the positions of coefficients to be tested.")
  }

  if (!is.numeric(q) & q[1] != "all") {
    stop("q should be a numeric vector indicating the m-quantiles for which
         coefficients are to be tested.")
  }

  if (q[1] != "all") {
    if (sum(q %in% object$q.values) != length(q)) {
      stop("You have chosen m-quantiles for which no coefficients have been estimated.
           Please estimate the coefficients before using mqwald(). See ?mquantreg.")
    }
  }

  # check if model has intercept

  intercept <- all(object$X[, 1] == 1)

  coef_model <- if (!is.matrix(object$coefficients)) {
                    matrix(object$coefficients[, -1],
                            nrow = length(object$q.values), byrow = TRUE)
                    } else {object$coefficients[, -1, drop = FALSE]}

  mq_model   <- if (!is.matrix(object$coefficients)) {
                  object$coefficients[1]
                  } else{object$coefficients[, 1]}

  # testing linear hypotheses, joint testing of coefficients
  if (isFALSE(each)) {
  p <- if (id.beta[1] == "all") {
    dim(coef_model)[2] - sum(intercept)
    } else {
      length(id.beta)
    }

  id.beta <- if (id.beta[1] == "all" & intercept == TRUE) {
    seq(from = 2, to = p + 1)
  } else if (id.beta[1] == "all" & intercept == FALSE) {
    seq(from = 1, to = p)
  } else {id.beta}

  mq <- if (q[1] == "all") {
    1:length(mq_model)
  } else {
    which(q %in% mq_model)
  }

  # create R matrix, if is.null in argument of function
  if (is.null(R)) {
    Rncol <- ncol(coef_model)
    Rnrow <- p
    R_test     <- matrix(rep(0, Rncol * Rnrow), nrow = Rnrow)
    test_tmp <- id.beta
    for (i in 1:nrow(R_test)) {
      R_test[i, test_tmp[i]] <- 1
    }
  } else {
    R_test <- R
  }


  if (is.null(H0)) {
    rrow <- p
    r    <- rep(0, times = rrow)
  } else {
    r <- H0
  }



  wald.test <- vector(mode = "numeric", length = length(mq))
  p.values  <- vector(mode = "numeric", length = length(mq))

  if (p > 1) {

    for (i in mq) {
    pos      <- which(i == mq)
    beta.q   <- coef_model[i, ]
    var.beta <- object$var.beta.matrix[ , , i]

    wald.test[pos] <- t(R_test %*% beta.q - r) %*% solve(R_test %*% var.beta %*% t(R_test)) %*% (R_test %*% beta.q - r)

    p.values[pos]  <- 1 - pchisq(wald.test[pos],p)

    }
  }

  if (p == 1) {
    for (i in mq) {
      pos      <- which(i == mq)
      beta.q   <- coef_model[i, id.beta]
      var.beta <- as.numeric(object$var.beta.matrix[id.beta, id.beta, i])

      wald.test[pos] <- (beta.q - r) * (1 / var.beta) * (beta.q - r)

      p.values[pos]  <- 1 - pchisq(wald.test[pos],p)
    }
  }
  res        <- list()
  class(res) <- "mqwald"

  res$call         <- cl
  res$coefficients <- object$coefficients
  res$varCov       <- object$var.beta.matrix
  res$df           <- rep(p, length(wald.test))
  res$statistic    <- wald.test
  res$p.values     <- p.values
  res$R            <- R_test
  res$r            <- r
  res$id.beta      <- id.beta
  res$q            <- q
  res$print_id     <- print_id


  colnames(res$R) <- colnames(object$coefficients)[-1]

  return(res)
  }

  # testing coefficients individually
  if(isTRUE(each)) {
    if(id.beta[1] == "all") {
      stop("Please choose which coefficients should be tested individually.
           The argument id.beta being set to 'all' is used for omnibus-testing the coefficients.")
    }

    if (!is.null(R)) {
      stop("No R matrix is needed. Coefficients are to be tested individually. Please set R to NULL.")
    }

    mq <- if (q[1] == "all") {
      1:length(mq_model)
    } else {
      which(mq_model == q)
    }

    if (is.null(H0)) {
      rrow <- length(id.beta)
      r    <- rep(0, times = rrow)
    } else {
      r <- H0
    }

    if (length(r) != length(id.beta)){
      stop("The length of the Null hypothesis vector and the number of coefficients to be tested don't coincide.
           Please make sure if you choose to test against values different from 0, that the Null hypothesis vector
           has the appropriate length.")
    }

    wald.test <- matrix(nrow = length(mq), ncol = length(id.beta))
    p.values  <- matrix(nrow = length(mq), ncol = length(id.beta))

    for (i in mq) {
      # i <- 1
      pos_row      <- which(i == mq)
      for (j in id.beta) {
        # j <- 2
        pos_col  <- which(j == id.beta)
        beta.q   <- coef_model[i, j]
        var.beta <- as.numeric(object$var.beta.matrix[j, j, i])

        wald.test[pos_row, pos_col] <- (beta.q - r[pos_col]) * (1 / var.beta) * (beta.q - r[pos_col])
        p.values[pos_row, pos_col]  <- 1 - pchisq(wald.test[pos_row, pos_col], 1)
        }
      }
    res        <- list()
    class(res) <- "mqwald"

    res$call         <- cl
    res$coefficients <- object$coefficients
    res$varCov       <- object$var.beta.matrix
    res$df           <- rep(1, length(wald.test))
    res$statistic    <- wald.test
    res$p.values     <- p.values

    res$R            <- NULL
    res$r            <- r
    res$id.beta      <- id.beta
    res$q            <- q
    res$print_id     <- print_id

    colnames(res$statistic) <- colnames(coef_model)[id.beta]
    # rownames(res$statistic) <- print_id$q
    colnames(res$p.values)  <- colnames(coef_model)[id.beta]

    return(res)

  }
}


#' Print Wald tests for M-Quantile regression coefficients
#'
#' The Print-method returns the tests for the chosen coefficients and m-quantiles in one table.
#'
#'@param x a mqwald object, produced by calling the mqwald-function.
#'@examples
#'df = simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.05)
#'fit  = mquantreg(formula = Y ~ x1 + x2,data = df, q  = 0.5, method = "continuous")
#'test = mqwald(fit)
#'print(test)
#'@export
print.mqwald <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)

  if (obj$print_id$each == FALSE & obj$print_id$id.beta[1] == "all" & is.null(obj$print_id$R) & is.null(obj$print_id$H0)) {
    qval <- if(obj$q[1] == "all") {
      obj$coefficients[, 1]
    } else {
      obj$q
    }
    cat("\nOmnibus test of significance for q-value(s)", paste(qval, collapse = ", "))
    cat("\nR:")
    print(obj$R)
    cat("\nH0:")
    print(obj$r)

  }

  if (obj$print_id$each == FALSE & (!is.null(obj$print_id$R) | !is.null(obj$print_id$H0))) {
    qval <- if(obj$q == "all") {
      obj$coefficients[, 1]
    } else {
      obj$q
    }

    Rprint <- if(is.null(obj$print_id$R)) {
      obj$R
    } else {
      obj$print_id$R
    }


    rprint <- if(is.null(obj$print_id$H0)) {
      obj$r
    } else {
      obj$print_id$H0
    }
    cat("\nTest of user specified hypotheses for q-value(s)", paste(qval, collapse = ", "))
    cat("\nR:")
    print(Rprint)
    cat("\nH0:")
    print(rprint)

  }

  if (obj$print_id$each == TRUE) {
    qval <- if(obj$q == "all") {
      obj$coefficients[, 1]
    } else {
      obj$q
    }

    rprint <- if(is.null(obj$print_id$H0)) {
      obj$r
    } else {
      obj$print_id$H0
    }
    cat("\nWald tests for single coefficients at q-value(s)", paste(qval, collapse = ", "))
    cat("\nH0:")
    print(rprint)

  }


  cat("\ntest-statistic:\n")
  print(obj$statistic)
  cat("\ndegrees of freedom\n")
  print(obj$df)
  cat("\np-values\n")
  print(obj$p.values)

}


#' Likelihood-Ratio test for M-quantile regression
#'
#' This function allows to test the Null hypothesis of regression coefficients
#' of an estimated m-quantile model to be jointly zero. This is accomplished by comparing
#' an estimated 'full' model with a user specified 'reduced' model, i.e. models have to be nested.
#'
#' @param object an object of class "mquantfit", usually, a result of a call to \link{mquantreg}
#' @param id.beta identifier for the coefficient(s) to be jointly tested. Expects vector of positions
#' in coefficients of \code{object}. Defaults to "all", which does not include the intercept.
#' @param q identifies the m-quantile(s) for which coefficients should be tested.
#' Expects vector of m-quantiles e.g. \code{c(0.3, 0.5, 0.7)}.Defaults to "all".
#' @param alpha significance level can be chosen to get a statement about significance, when printing.
#'
#' @return the function mqLRT returns the the result of
#'
#' @param coefficents the estimated coeffients of the full model to be tested.
#' @param statistic the test statistic
#' @param p.values the p-values for the wald tests.
#' @param df the degrees of freedom of the LR-tests.
#'
#'
#'
#' @examples
#' df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                   measurement.error = 0.05)
#'
#' fit <- mquantreg(formula = Y ~ x1 + x2, data = df, q  = c(0.5, 0.7), method = "continuous")
#'
#' hyp_test <- mqLRT(fit, id.beta = "all", q = "all")
#' print(hyp_test)
#' @export
mqLRT=function(object, id.beta = "all", q = "all", alpha = 0.05){
  cl <- match.call()
  # object <- fit
  # id.beta = c(1:3)
  # q = "all"

  if(object$method != "continuous"){
    stop("Likelihood-Ratio hypothesis testing is only implemented for continuous dependent variables
         estimated with the likelihood method.")
  }

  if(!is.numeric(id.beta) & id.beta[1] != "all"){
    stop("id.beta must be an vector indicating the positions of coefficients to be tested.")
  }

  if(!is.numeric(q) & q[1] != "all"){
    stop("q should be a numeric vector indicating the m-quantiles for which
         coefficients are to be tested.")
  }

  if(q[1] != "all") {
    if(sum(q %in% object$q.values) != length(q)) {
      stop("You have chosen m-quantiles for which no coefficients have been estimated.
           Please estimate the coefficients before using mqwald(). See ?mquantreg.")
    }
  }

  intercept <- all(object$X[, 1] == 1)

  if(intercept == FALSE & id.beta[1] == "all") {
    stop("The model to be tested has no intercept and 'all' coefficients were chosen
         to be tested. Please either choose less coefficients to be tested or
         reestimate the model with an intercept before setting id.beta to 'all'.")
  }

  coef_model <- if(!is.matrix(object$coefficients)) {
    matrix(object$coefficients[, -1],
           nrow = length(object$q.values), byrow = TRUE)
  } else {object$coefficients[, -1, drop = FALSE]}

  mq_model   <- if(!is.matrix(object$coefficients)){
    object$coefficients[1]
  } else{object$coefficients[, 1]}

  # number of estimated coefficients full model and restricted model
  p_full <- ncol(coef_model)
  p_red  <- ifelse(id.beta[1] == "all", 1, p_full - length(id.beta))

  n <- nrow(residuals(object))

  # design matrix restricted model
  x1 <- if(id.beta[1] == "all") {
    as.matrix(object$X[, 1])
    } else {
    object$X[, - c(id.beta), drop = FALSE]
    }

  if(id.beta[1] == "all") colnames(x1) <- "(Intercept)"
  if(any(colnames(x1) %in% "(Intercept)")) colnames(x1)[which(colnames(x1) == "(Intercept)")] <- "1"



  q      <- if(q[1] == "all"){sort(object$q.values)} else {sort(q)}
  LRtest <- vector(mode = "numeric", length = length(q))
  prob   <- vector(mode = "numeric", length = length(q))

  for (i in 1:length(q)){

    # model estimation full and reduced
    if(class(object$formula) == "formula") object$formula <- paste(as.character(object$formula)[c(2, 1, 3)], collapse = " ")

    tmp_full <- mquantreg(formula = object$formula, data = object$data, q = q[i],
                          method = "continuous", k = object$k)
    formula_tmp <- paste0(strsplit(object$formula, " ")[[1]][1]," ~ ", paste(colnames(x1), collapse = " + "))
    if(colnames(x1)[1] == "(intercept)" & length(colnames(x1)) == 1) formula_tmp <- paste0(strsplit(object$formula, " ")[[1]][1]," ~ 1")
    dat_tmp     <- as.data.frame(cbind(object$Y, x1))
    colnames(dat_tmp)[1] <- strsplit(object$formula, " ")[[1]][1]
    tmp_red  <- mquantreg(formula = formula_tmp, data = dat_tmp, q = q[i],
                          method = "continuous", k = object$k)

    # residuals full and reduced
    res_full <- tmp_full$resid
    res_red  <- tmp_red$resid

    # full model
    sm <- tmp_full$scale
    w  <- (object$k * (abs(res_full / sm) - object$k / 2)) * ((res_full / sm) <=
                - object$k | (res_full / sm) >= object$k) + 0.5 * (res_full / sm)^2 *
                (-object$k < (res_full / sm) & (res_full / sm) < object$k)

    ww              <- 2 * (1 - q[i]) * w
    ww[res_full> 0] <- 2 * q[i] * w[res_full > 0]
    w               <- ww
    V_full          <- sum(w)

    # reduced model
    sm2 <- sm
    object$k2 <- object$k
    w2 <- (object$k2 * (abs(res_red / sm2) - object$k2 / 2)) * ((res_red / sm2) <=
                - object$k2 | (res_red / sm2) >= object$k2) + 0.5* (res_red / sm2)^2 *
                (-object$k2 < (res_red / sm2) & (res_red / sm2) < object$k2)

    ww2             <- 2 * (1 - q[i]) * w2
    ww2[res_red> 0] <- 2 * q[i] * w2[res_red > 0]
    w2              <- ww2
    V_red           <- sum(w2)

    Epsi   <- sum((2 * (q[i] * (0 <= res_full / sm & res_full / sm <= object$k) +
                (1 - q[i]) * (- object$k <= res_full / sm & res_full / sm<0)))) / (n)
    res_q                <- psi.huber(res_full / sm, k = object$k) * (res_full / sm)
    res_qw               <- 2 * (1 - q[i]) * res_q
    res_qw[res_full > 0] <- 2 * q[i] * res_q[res_full > 0]
    res_q                <- res_qw

    Epsi2 <- (sum((res_q)^2) / (n - p_red))

    # estimation of test statistic and p-value
    df        <- p_full - p_red
    LRtest[i] <- -2 * (Epsi / (Epsi2)) * (V_full - V_red)
    prob[i]   <- 1 - pchisq(LRtest[i], df)


  }

  beta <- if(id.beta[1] == "all") {seq(3, p_full + 1)} else {id.beta + 1}
  res        <- list()
  class(res) <- "mqLRT"

  res$call         <- cl
  res$coefficients <- object$coefficients[object$q.values %in% q, c(1,beta)]
  res$statistic    <- LRtest
  res$p.values     <- prob
  res$df           <- df
  res$alpha        <- alpha
  return(res)
}

#' Print Likelihood-Ratio test results for M-Quantile regression coefficients
#'
#' The Print-method returns the tests for the chosen coefficients and m-quantiles in one table.
#'
#'@param x a mqLRT object, produced by calling the mqLRT-function.
#'@examples
#'df = simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.05)
#'fit  = mquantreg(formula = Y ~ x1 + x2, data = df, q  = 0.5, method = "continuous")
#'test = mqLRT(fit)
#'print(test)
#'@export
print.mqLRT <- function(x, verbose = TRUE, ...){
  # x <- test
  # verbose = TRUE
  cat("Call:\n")
  print(x$call)
  coef <- x$coefficients
  cat("\nEstimated coefficients to be restricted to zero in reduced model:\n")
  print(coef, digits = 4)

  coef_res <- if(!is.matrix(coef)){coef[1]} else {coef[, 1]}
  result <- matrix(c(coef_res, x$statistic, rep(x$df, length(x$p.values)), x$p.values), ncol = 4, byrow = FALSE)
  colnames(result) <- c("q-value(s)", "Chi^2", "df", "p-value(s)")

  cat("\nResult of Likelihood-Ratio test(s):\n")
  print(result)

  if(verbose == TRUE){
    sig <- x$p.values <= x$alpha
    sig <- if(!is.matrix(x$coefficients)){x$coefficients[1][sig]} else {x$coefficients[, 1][sig]}
    alpha <- x$alpha
    if(length(sig) == 0){
      cat("\nThere were no significant Likelihood-Ratio tests at the chosen alpha level of ", alpha, ".\n")
    } else {
    cat("\nThe likelihood-ratio test gives significant results at a level alpha of",
        alpha, "for the m-quantile(s)", paste(sig, collapse = ", "),".")
    }
  }
}

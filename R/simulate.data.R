#' Simulate Data
#'
#' The \code{simulate_data}-function generates data, which can be used for testing and exploring m-quantile regression.
#' The data can be contaminated with outliers for testing the robust properties of m-quantile regression.
#'
#'@param n numeric, number of observations to generate.
#'@param real.betas vector with true beta used for simulation.
#'@param X dataframe with x-values, if left empty x-values will be generated as well.
#'@param response.type string indicating what type of error term the data should have.
#' Following distribution are available:
#'\itemize{
#'\item \code{"continuous.normal"} - normally distributed data with symmetric outlier contamination.
#'\item \code{"continuous.normal.posout"} - normally distributed data with only positive outlier contamination.
#'\item \code{"continuous.lognormal"} - log-normally distributed data.
#'\item \code{"count.nb"} - negative binomially distributed data.
#'\item \code{"count.poisson"} -  poisson distributed data.
#'\item \code{"binary"} - binary data
#'}
#'@param measurement.error numeric, controlling what proportion of the data is contaminated with outliers.
#'@param sd numeric setting the standard deviation of the error term.
#'@param error.strength numeric controlling the strangth of outliers.
#'@return dataframe
#'@examples
#'data_sim <- simulate_data(300,
#'                         real.betas = c(0, 5),
#'                         X = FALSE,
#'                         response.type = "continuous.normal",
#'                         measurement.error = 0.05,
#'                         sd = 1,
#'                         error.strength = 3)
#'@export


simulate_data = function(n,
                         real.betas = c(0, 1),
                         X = FALSE,
                         response.type = "continuous.normal",
                         measurement.error = 0.05,
                         sd = 1,
                         error.strength = 3)
{
  ##############################################
  ## Check if n is valid
  if (n < 0)
    stop("n must be a number larger than 0")

  ##############################################
  ## Check if valid method is selected
  if (!(response.type %in% c("continuous.normal", "continuous.lognormal", "continuous.normal.posout", "count.nb", "count.poisson", "binary")))
  {
    stop(
      "Chosen response type is not implemented. Avaiable methods are: continuous.normal, continuous.lognormal, binary,  count.nb, count.poisson, continuous.normal.posout"
    )
  }

  ##############################################
  ## Check if valid m-quantile q is choosen
  if (any(measurement.error < 0) | any (measurement.error > 1))
    stop(
      "Choosen proportion of measurement error is not valid. q must be greater that 0 and lesser than 1."
    )

  ##############################################
  ## Check if X is passed
  ind.number = length(real.betas) - 1
  if (any(X == FALSE)) {
  if(ind.number == 0){
      X = rep(1, n)
      dim(X) = c(n, 1)
      X = data.frame(X)
      colnames(X) = c("intercept")
    }
    else{
      X = c(rep(1, n), rnorm(n * ind.number, mean = 0, sd = 1))
      dim(X) = c(n, length(real.betas))
      X = data.frame(X)
      colnames(X) = c("intercept", paste0("x", 1:ind.number))
    }

  }
  data = as.data.frame(X)

  ##############################################
  ## Check if measurement error is valid
  if (any(measurement.error < 0) | any (measurement.error > 1))
    stop(
      "Chosen proportion of measurement error is not valid. Proportion must be greater that 0 and lesser than 1."
    )

  if (response.type == "continuous.normal") {
    data$Y = rnorm(n, mean = as.matrix(X) %*% real.betas, sd = sd)

    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    n.clean = n - n.error
    error = sd(data$Y) * error.strength
    error.var = sample(c(rep(error, n.error), rep(0, n.clean)))
    data$Y = data$Y +  replicate(n, sample(c(-1, 1), 1)) * error.var
  }

  if (response.type == "continuous.lognormal") {
    # Generate Y from negative binomial distribution
    errors <- rlnorm(n, meanlog=0, sdlog= sd)
    errors <- errors - 1.65  # this centers the error distribution on 0
    data$Y = as.matrix(X) %*% real.betas + errors


    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    n.clean = n - n.error
    error = sd(data$Y) * error.strength
    error.var = sample(c(rep(error, n.error), rep(0, n.clean)))
    data$Y = data$Y +  replicate(n, sample(c(-1, 1), 1)) * error.var
  }

  if (response.type == "continuous.normal.posout") {
    if (any(measurement.error >=0.5))
      stop(
        "Choosen proportion of measurement error is not valid for simualation with positive outliers. The proportion of measurement errors must be below 0.5"
      )

    eta = as.matrix(X) %*% real.betas
    errorterm = rnorm(n, mean = 0, sd = sd)
    data$Y = eta + errorterm

    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    poserror = errorterm > 0
    n.clean = poserror - n.error
    n.poserror = sum(poserror)
    error = sd(data$Y) * error.strength
    error.var = rep(0, n)

    error.var[poserror] = sample(c(rep(error, n.error), rep(0, n.poserror - n.error)))
    data$Y = data$Y +  error.var
  }

  if (response.type == "count.poisson") {
    # Generate Y with normal distributed errorterm
    eta = rnorm(n, mean = as.matrix(X) %*% real.betas, sd = sd)
    data$Y = rpois(n, lambda = exp(eta))

    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    n.clean = n - n.error
    error = floor(sd(data$Y) * error.strength)
    error.var = sample(c(rep(error, n.error), rep(0, n.clean)))
    data$Y = data$Y + error.var
  }

  if (response.type == "count.nb") {
    # Generate Y with normal distributed errorterm
    mu = exp(rnorm(n, mean = as.matrix(X) %*% real.betas, sd = sd))
    theta = 4
    data$Y = rnegbin(n, mu = mu, theta = theta)

    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    n.clean = n - n.error
    error = floor(sd(data$Y) * error.strength)
    error.var = sample(c(rep(error, n.error), rep(0, n.clean)))
    data$Y = data$Y +  error.var
  }

  if (response.type == "binary") {
    # Generate Y with normal distributed errorterm
    eta = rnorm(n, mean = as.matrix(X) %*% real.betas, sd = sd)
    p = exp(eta) * (1 + exp(eta))^(-1)
    data$Y = rbinom(n, size = 1, prob = p)#Bernoulli distribution

    # Add measurement error to response variable with
    n.error = floor(n * measurement.error)
    n.clean = n - n.error
    logit_fit = glm(formula = "Y ~ .", data = data, family=binomial(link="logit"))
    props = logit_fit$fitted.values
    percentile = ecdf(props)
    outlier = percentile(props) > 1 - measurement.error
    data$Y[outlier] = 1
  }


  return(data)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyp_testing.R
\name{print.mqwald}
\alias{print.mqwald}
\title{Print Wald tests for M-Quantile regression coefficients}
\usage{
\method{print}{mqwald}(obj, ...)
}
\arguments{
\item{x}{a mqwald object, produced by calling the mqwald-function.}
}
\description{
The Print-method returns the tests for the chosen coefficients and m-quantiles in one table.
}
\examples{
df = simulate_data(n = 1000,
                  real.betas = c(0.1, 0.3, 0.1 ),
                  response.type = "continuous.normal",
                   measurement.error = 0.05)
fit  = mquantreg(formula = Y ~ x1 + x2,data = df, q  = 0.5, method = "continuous")
test = mqwald(fit)
print(test)
}

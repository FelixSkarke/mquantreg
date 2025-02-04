% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{ggplot.prop.huber.q}
\alias{ggplot.prop.huber.q}
\title{Proportion of residuals huberised over q}
\usage{
\method{ggplot}{prop.huber.q}(fit)
}
\arguments{
\item{fit}{a mquantfit object, produced by calling the mquantreg-function.}
}
\description{
This plot shows the effect of k on the amount of observation being classified as outliers and therefore downweighted.
}
\examples{
library(mq1)

df <- simulate_data( n = 100,
                  real.betas = c(0.1, 0.3, 0.1 ),
                  response.type = "continuous.normal",
                   measurement.error = 0.5)

fit <- mquantreg(formula = "Y ~ x1 + x2",data = df, q  = 0.5, method = "continuous")
plot(fit, plottype = "prop.huber.q")
}

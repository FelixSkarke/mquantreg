#' Backward stepwise variable selection using qscores.
#'
#' The function implements a backward variable selection method based on the correlation coefficient
#' of the models respective qscores.

#'
#'@param object object containing regression model fitted with mquantreg
#'@param scope 	defines the range of models examined in the stepwise search. This should be either a single formula,
#'or a list containing the component lower, a formula. See the details for how to specify the formula and how the argument is used.
#'@param threshold number strictly between 0 and 1. variables get dropped, if correlation of qscores
#'with model from step before exceeds the threshold.
#'@param trace if positive, information is printed during the running of stepmq.
#'@param keep a filter function whose input is a fitted model object and the associated qscores, and whose output is arbitrary.
#'Typically keep will select a subset of the components of the object and return them. The default is not to keep anything.
#'@param steps the maximum number of steps to be considered. The default is 1000 (essentially as many as required). It is
#'typically used to stop the process early.

#'@return See \code{\link{mquantreg}} for details.
#'@references
#'Dawber, J. (2017). Advances in M-quantile estimation, University of Wollongong Thesis Collections
#'@examples
#'library(mq1)
#'df <- simulate_data(n = 1000,
#'                   real.betas = c(0.1, 0.3, 0.1 ),
#'                   response.type = "continuous.normal",
#'                    measurement.error = 0.5)
#'fit <- mquantreg(formula = Y ~ x1 + x2, data = df, q  = 0.5, method = "continuous")
#'stepmq(fit)
#'@export
stepmq <- function(object, scope, threshold = 0.95, trace = 1, keep = NULL, steps = 1000) {
  # browser()
   model_sat        <- object
   model_sat$call$q <- 0.5
   thresh           <- threshold


   cut.string <- function(string) {
     if (length(string) > 1L)
       string[-1L] <- paste0("\n", string[-1L])
     string
   }

   re.arrange <- function(keep) {
     namr <- names(k1 <- keep[[1L]])
     namc <- names(keep)
     nc <- length(keep)
     nr <- length(k1)
     array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
                                                            namc))
   }

   step.results <- function(models, fit, object) {
     change  <- sapply(models[seq(2)], "[[", "change")
     rdf     <- sapply(models[seq(2)], "[[", "df.resid")
     ddf     <- c(NA, diff(rdf))
     qsc     <- sapply(models[seq(2)], "[[", "qscores")
     heading <- c("Stepwise Model Path \nAnalysis of qscores correlations",
                  "\nInitial Model:", deparse(formula(object)),
                  "\nFinal Model:", deparse(formula(fit)), "\n")
     aoqs    <- data.frame(Step = I(change), Df = ddf,
                           `Resid. Df` = rdf,
                           check.names = FALSE)

     attr(aoqs, "heading") <- heading
     fit$qscores_aoqs      <- aoqs
     fit
   }

   Terms <- terms(object)
   object$call$formula <- object$formula <- Terms


   if (missing(scope)) {
     fdrop <- numeric()
     # fadd <- attr(Terms, "factors")
   } else {
     if (is.list(scope)) {
       fdrop <- if (!is.null(fdrop <- scope$lower))
         attr(terms(update.formula(object, fdrop)), "factors")
       else numeric()
     }
   }

   models <- vector("list", steps)

   if (!is.null(keep)) keep.list <- vector("list", steps)

   n            <- nobs(object, use.fallback = TRUE)
   fit          <- object
   qscores_init <- qscores(fit)
   df           <- n - length(coef(fit)[-1L])
   nm           <- 1

   if (trace) {
     cat("Starting model:", "\n",
         cut.string(deparse(formula(fit))), "\n\n",
         sep = "")
     flush.console()
   }
   models[[nm]] <- list(qscores = qscores_init, df.resid = df, change = "") # neu: change = 0

   if (!is.null(keep))
     keep.list[[nm]] <- keep(fit, qscores_init)

   while (steps > 0) {
     steps        <- steps - 1
     qscores_step <- qscores_init
     ffac         <- attr(Terms, "factors")
     scope        <- factor.scope(ffac, list(drop = fdrop))
     aoqs         <- NULL
     change       <- NULL

     if (length(scope$drop)) {
       aoqs <- drop1.mquantreg(object = fit, scope = scope$drop, trace = trace)
       rn   <- row.names(aoqs)

       row.names(aoqs) <- c(rn[1L], paste("-", rn[-1L]))

       if (any(aoqs$Df == 0, na.rm = TRUE)) {
         zdf <- aoqs$Df == 0 & !is.na(aoqs$Df)
         change <- rev(rownames(aoqs)[zdf])[1L]
       }
     }
     if (is.null(change)) {
       attr(aoqs, "heading") <- NULL

       nzdf <- if (!is.null(aoqs$Df)) {
         aoqs$Df != 0 | is.na(aoqs$Df)
       }

       aoqs <- aoqs[nzdf, ]

       if (is.null(aoqs) || ncol(aoqs) == 0)
         break

       nc <- match("Corr.", names(aoqs))
       nc <- nc[!is.na(nc)][1L]
       o <- order(aoqs[, nc], decreasing = TRUE)

       if (trace)
         print(aoqs[o, ])

       if (nrow(aoqs) == 1)
         break

       if (aoqs[o[2L], nc] < threshold)
         break

       change <- rownames(aoqs)[o[2L]]
     }

     fit <- update(fit, paste("~ .", change), evaluate = FALSE)
     fit <- eval.parent(fit)
     nnew <- nobs(fit, use.fallback = TRUE)
     if (all(is.finite(c(n, nnew))) && nnew != n)
       stop("number of rows in use has changed: remove missing values?")

     Terms <- terms(fit)

     nm <- nm + 1
     models[[nm]] <- list(qscores = qscores(fit), df.resid = n -
                            length(coef(fit)[-1]), change = change)
     if (!is.null(keep))
       keep.list[[nm]] <- keep(fit, qscores(fit))
   }

   if (!is.null(keep))
     fit$keep <- re.arrange(keep.list[seq(nm)])

   step.results(models = models[seq(nm)], fit, object)


 }


drop1.mquantreg <- function (object, scope, trace = FALSE) {
  tl <- attr(terms(object), "term.labels")
  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>",
            scope), c("df", "corr.")))

  n0 <- nobs(object, use.fallback = TRUE)
  df0 <- n0 - (dim(coef(object))[2] - 1)
  env <- environment(formula(object))
  qscores_nfit <- vector("list", ns + 1)
  qscores_nfit[[1]] <- qscores(object)

  for (i in seq_len(ns)) {
    tt <- scope[i]
    if (trace > 1) {
      cat("trying -", tt, "\n", sep = "")
      flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . -",
                                            tt)), evaluate = FALSE)
    nfit <- eval(nfit, envir = env)

    qscores_nfit[[i + 1]] <- qscores(nfit)
    nnew                  <- nobs(nfit, use.fallback = TRUE)

    if (all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")

    dfnew         <- nnew - (dim(coef(nfit))[2] - 1)
    ans[i + 1, 1] <- dfnew
  }

  for( i in seq_along(qscores_nfit)){
    ans[i, 2] <- cor(qscores_nfit[[1]], qscores_nfit[[i]])
  }

  ans[1L, 1L] <- df0
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  ans[, 1L] <- dfs
  aoqs <- data.frame(Df = ans[,1], Corr. = ans[, 2])

  head        <- c("Single term deletions", "\nModel:",
                   deparse(formula(object)))
  class(aoqs) <- c("corr.qscore", "data.frame")
  attr(aoqs, "heading") <- head
  aoqs
}




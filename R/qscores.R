main.qscores <- function(y, yhatq, qvals) {
  n <- length(y)

  if (nrow(yhatq) != n){
    stop("y-dimensions do not agree")
  }

  Q <- length(qvals)

  if (ncol(yhatq) != Q){
    stop("q-dimensions do not agree")
  }

  qvec   <- c(qvals[1] / n, qvals, (n - 1 + qvals[Q]) / n)
  qscore <- rep(0, n)

  set.seed(123)
  for (i in 1:n) {
    nq <- (1:Q)[abs(y[i] - yhatq[i, ]) == min(abs(y[i] - yhatq[i, ]))]
    lq <- length(nq)

    if (lq > 1) {
      if (y[i] > yhatq[i, nq[lq]]) {
        nq <- nq[lq]
      } else {
        nq <- nq[1]
      }

    }

    q1 <- qvec[nq + 1]

    if (y[i] > yhatq[i, nq]) {
      q2 <- qvec[nq + 2]
      if (nq == Q) {
        fac = runif(1)
      } else {
        fac = (y[i] - yhatq[i, nq]) / (yhatq[i, nq + 1] - yhatq[i, nq])
      }
      qscore[i] = q1 * (1 - fac) + q2 * fac

    } else {
      q2 <- qvec[nq]
      if (nq == 1) {
        fac = runif(1)
      } else{
        fac = (y[i] - yhatq[i, nq - 1]) / (yhatq[i, nq] - yhatq[i, nq - 1])
      }
      qscore[i] = q1 * fac + q2 * (1 - fac)
    }
  }
  qscore[qscore > 1] <- 1
  qscore[qscore < 0] <- 0
  qscore
}

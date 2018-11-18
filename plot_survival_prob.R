plot_survival_prob<-function (x, from = NULL, to = NULL, range = NULL, covariates = "mean", 
          legend.pos = NULL, xlab = "Time", ylab = "Fitted survival probability", 
          lwd = 1, ylim = NULL, ...) 
{
  if (!inherits(x, "msm")) 
    stop("expected x to be a msm model")
  if (is.null(from)) 
    from <- transient.msm(x)
  else {
    if (!is.numeric(from)) 
      stop("from must be numeric")
    if (any(!(from %in% 1:x$qmodel$nstates))) 
      stop("from must be a vector of states in 1, ..., ", 
           x$qmodel$nstates)
  }
  if (is.null(to)) 
    to <- max(absorbing.msm(x))
  else {
    if (!is.numeric(to)) 
      stop("to must be numeric")
    if (!(to %in% absorbing.msm(x))) 
      stop("to must be an absorbing state")
  }
  if (is.null(range)) 
    rg <- range(model.extract(x$data$mf, "time"))
  else {
    if (!is.numeric(range) || length(range) != 2) 
      stop("range must be a numeric vector of two elements")
    rg <- range
  }
  if (is.null(ylim)) 
    ylim <- c(0,1)
  else {
    ylim <- ylim
  }
  timediff <- (rg[2] - rg[1])/50
  times <- seq(rg[1], rg[2], timediff)
  pr <- numeric()
  cols <- rainbow(length(from))
  for (t in times) pr <- c(pr, pmatrix.msm(x, t, times[1], 
                                           covariates)[from[1], to])
  plot(times, 1 - pr, type = "l", xlab = xlab, ylab = ylab, 
       lwd = lwd, ylim = ylim, lty = 1, col = cols[1], ...)
  lt <- 2
  for (st in from[-1]) {
    pr <- numeric()
    for (t in times) pr <- c(pr, pmatrix.msm(x, t, times[1], 
                                             covariates)[st, to])
    lines(times, 1 - pr, type = "l", lty = lt, lwd = lwd, 
          col = cols[lt], ...)
    lt <- lt + 1
  }
  if (!is.numeric(legend.pos) || length(legend.pos) != 2) 
    legend.pos <- c(max(times) - 15 * timediff, 1)
  legend(legend.pos[1], legend.pos[2], legend = paste("From state",from), lty = seq(lt - 1), col = cols, lwd = lwd)
  invisible()
}
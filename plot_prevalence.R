plot_prevalence<-function (x, mintime = NULL, maxtime = NULL, timezero = NULL, 
          initstates = NULL, interp = c("start", "midpoint"), censtime = Inf, 
          subset = NULL, covariates = "population", misccovariates = "mean", 
          piecewise.times = NULL, piecewise.covariates = NULL, xlab = "Year", 
          ylab = "Prevalence (%)", lwd.obs = 2, lwd.exp = 2, lty.obs = 1, 
          lty.exp = 2, col.obs = "blue", col.exp = "red", legend.pos = NULL, 
          startyr = NULL, ...) 
{
  if (!inherits(x, "msm")) 
    stop("expected x to be a msm model")
  time <- model.extract(x$data$mf, "time")
  if (is.null(mintime)) 
    mintime <- min(time)
  if (is.null(maxtime)) 
    maxtime <- max(time)
  t <- seq(mintime, maxtime, length = 100)
  obs <- msm:::observed.msm(x, t, interp, censtime, subset)
  expec <- msm:::expected.msm(x, t, timezero = timezero, initstates = initstates, 
                        covariates = covariates, misccovariates = misccovariates, 
                        piecewise.times = piecewise.times, piecewise.covariates = piecewise.covariates, 
                        risk = obs$risk, subset = subset, ci = "none")[[2]]
  states <- seq(length = x$qmodel$nstates)
  S <- length(states)
  ncols <- ceiling(sqrt(S))
  nrows <- if (floor(sqrt(S))^2 < S && S <= floor(sqrt(S)) * 
               ceiling(sqrt(S))) 
    floor(sqrt(S))
  else ceiling(sqrt(S))
  par(mfrow = c(nrows, ncols))
  for (i in states) {
    plot(startyr+t, obs$obsperc[, i], type = "l", ylim=c(max(0,min(obs$obsperc[, i])-30),min(100,max(obs$obsperc[, i])+10)),
         xlab = xlab, ylab = ylab, lwd = lwd.obs, lty = lty.obs, 
         col = col.obs, main = rownames(x$qmodel$qmatrix)[i], xaxt = "n", 
         ...)
    axis(1,at = c(2003,2004), labels = c("2003","2004"))
    lines(startyr+t, expec[, i], lwd = lwd.exp, lty = lty.exp, col = col.exp)
  }
  if (!is.numeric(legend.pos) || length(legend.pos) != 2) 
    legend.pos <- c(0.4 * maxtime, 40)
  legend(x = legend.pos[1], y = legend.pos[2], legend = c("Observed", 
                                                          "Expected"), lty = c(lty.obs, lty.exp), lwd = c(lwd.obs, 
                                                                                                          lwd.exp), col = c(col.obs, col.exp))
  invisible()
}
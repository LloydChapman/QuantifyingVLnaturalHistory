plot_obs_num_in_each_state<-function (x, mintime = NULL, maxtime = NULL, startyr = NULL,
          timezero = NULL, initstates = NULL, interp = c("start", "midpoint"), 
          censtime = Inf, subset = NULL, covariates = "population", 
          misccovariates = "mean", piecewise.times = NULL, piecewise.covariates = NULL, 
          xlab = "Time", ylab = "Number in state", lwd.obs = 2, lty.obs = 1, 
          legend.pos = NULL, ...) 
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
  states <- seq(length = x$qmodel$nstates)
  S <- length(states)
#   ncols <- 2
#   nrows <- 1 
#   par(mfrow = c(nrows, ncols))
  plot(0,type = "n", xlim=c(startyr+min(t),startyr+max(t)), ylim = c(0, 1000), 
       xlab = xlab, ylab = ylab, xaxt = "n")
  axis(1,at = c(2003,2004), labels = c("2003","2004"))
  # axis(1,at = c(1999,2000,2001,2002,2003,2004), labels = c("1999","2000","2001","2002","2003","2004"))
  col.obs<-c("green","yellow","red","blue","black")
  for (i in states) {
#     plot(startyr+t, obs$obstab[, i], type = "l", ylim = c(600, 1000), 
#          xlab = xlab, ylab = ylab, lwd = lwd.obs, lty = lty.obs, 
#          col = col.obs, main = rownames(x$qmodel$qmatrix)[i], 
#          ...)
    lines(startyr+t, obs$obstab[, i], lwd = lwd.obs, lty = lty.obs, 
          col = col.obs[i])
    # lines(t, expec[, i], lwd = lwd.exp, lty = lty.exp, col = col.exp)
  }
  if (!is.numeric(legend.pos) || length(legend.pos) != 2) 
    legend.pos <- c(startyr + 0.84 * maxtime, 600)
  # legend.pos <- c(startyr + 0.1 * maxtime, 600)
  legend(x = legend.pos[1], y = legend.pos[2], 
         legend = c("Susceptible","Asymptomatic","KA","Recovered/Dormant","Dead"), 
         lty = lty.obs, lwd = lwd.obs, col = col.obs)
  plot(0,type = "n", xlim=c(startyr+min(t),startyr+max(t)), ylim = c(0, 30), 
       xlab = xlab, ylab = ylab, xaxt = "n")
  axis(1,at = c(2003,2004), labels = c("2003","2004"))
  # axis(1,at = c(1999,2000,2001,2002,2003,2004), labels = c("1999","2000","2001","2002","2003","2004"))
  for (i in c(3,5)){
    lines(startyr+t, obs$obstab[, i], lwd = lwd.obs, lty = lty.obs, 
          col = col.obs[i])
  }
  legend(x = 2002.3, y = 30, 
         legend = c("KA","Dead"), 
         lty = lty.obs, lwd = lwd.obs, col = col.obs[c(3,5)])
  #   legend(x = 2001, y = 30, 
  #          legend = c("KA","Dead"), 
  #          lty = lty.obs, lwd = lwd.obs, col = col.obs[c(3,5)])
    invisible()
}
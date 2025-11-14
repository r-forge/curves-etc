#### called by sim.est() and  simestgcv()  in ./SimEst.r & ./SimEstGcv.r
##
## help() doc in ../man/smooth.pen.reg.Rd
smooth.pen.reg <- function (x, y, lambda, w = NULL, agcv = FALSE, agcv.iter = 100, ...)
{
	if (!all(is.finite(c(x, y))))
	  stop("missing or infinite values in inputs are not allowed")
	if ((n <- length(x)) != length(y))
	  stop("'x' and 'y' must have same length.")
	if (is.null(w))
            w <- rep_len(1, n)
	else {
            if(n != length(w)) stop("lengths of 'x' and 'w' must match")
            if(any(w <= 0))   stop("all weights should be positive")
	}
	if (!is.finite(lambda) || lambda < 0)
	  stop("'lambda' must be non-negative and finite")
        flag <- if(agcv) 1L else 0L
	## print(flag)
        if(is.unsorted(x)) { # sort along x
            io <- order(x)
            x <- x[io]
            y <- y[io]
            if(!is.null(w))
                w <- w[io]
        }
	h <- diff(x)
	Qty <- diff(diff(y)/h) # FIXME use central diffs (apart from boundary)
	out <- .C(spen_egcv, n = as.integer(n), x = as.double(x), y = as.double(y),
                  w = as.double(w),  h = as.double(h), Qty = as.double(Qty),
                  lm = as.double(c(lambda,0)), yhat = double(n), iter = as.integer(agcv.iter),
                  EGCV = flag, agcv = as.double(0))
        ## print(out$EGCV)
	yhat <- out$yhat
	spFn <- splinefun(x, yhat)# mainly to compute 1st deriv
	deriv <- unname(spFn(x, deriv = 1L))
	ret <- list(call = match.call(),
                    x.values = x, y.values = y, fit.values = yhat, w = w,
	            deriv = as.vector(deriv), residuals = y - yhat, iter = 1L,
	            convergence = 0,
                    minvalue = out$lm[2],
                    agcv.score = out$agcv,
                    splinefun = spFn)
	class(ret) <- "smooth.pen.reg"
        ret
}

print.smooth.pen.reg <- function(x, digits = getOption("digits"), ...) {
  cat("Call: \t"); print(x$call)
  cat(sprintf("Minimum Criterion Value Obtained: %s\n", format(x$minvalue, digits=digits)))
  cat(sprintf("Number of Iterations: %d\n", x$iter))
  cat(sprintf("Convergence Status:   %d\n", x$convergence))
  invisible(x)
}

plot.smooth.pen.reg <- function(x, diagnostics = TRUE,
                                ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
                                pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...)
{
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  if(diagnostics) {
      plot.window(c(0,7), c(0,7))
      par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0)) -> op; on.exit(par(op))
      plot(xx,yx, xlab = 'x', ylab=ylab, pch=pch, cex=cex,
           main = "Smoothing Spline using\n Penalized Least Squares", ...)
      lines(xx, fitx, lwd=lwd, col = col2)
      plot(fitx,resx,xlab = 'Fitted Values',ylab = "Residuals",pch=pch, main = "Fitted vs Residuals", ...)
      abline(h = 0.0, lty = ablty, col = col2)
      plot(yx,fitx,xlab = "Actual Values",ylab = "Fitted Values",pch=pch, main = "Actual vs Fitted", ...)
      abline(a = 0, b = 1, lty = ablty, col = col2)
      qqnorm(resx)
      qqline(resx)
  } else {
      plot.window(c(0,7), c(0,7))
      par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.3,.5,0)) -> op; on.exit(par(op))
      plot(xx,yx, xlab = 'x', ylab=ylab, pch=pch, cex=cex,
           main = "Smoothing Spline using\n Penalized Least Squares", ...)
      lines(xx, fitx, lwd=lwd)
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.smooth.pen.reg <- function(object, newdata = NULL, deriv = 0, ...){
    stopifnot("deriv must either be 0 or 1!" = (deriv == 0 || deriv == 1))
    if(!is.null(newdata))
        object$splinefun(as.vector(newdata), deriv = deriv)
    else {
        warning("No 'newdata' found and so using input 'x' values")
        if(deriv == 0)
            object$fit.values
        else # (deriv == 1)
            object$deriv
    }
}

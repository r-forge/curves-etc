cvx.pen.reg <- function(x, y, lambda, w = NULL, tol = 1e-05, maxit = 1000)
{
	if(length(x) != length(y))
		stop("'x' and 'y' must have same length.")
	if (!all(is.finite(c(x, y))))
		stop("missing or infinite values in inputs are not allowed")
	n <- length(x)
        if(n < 5)
            stop("Number of samples must be at least 5.")
	if (has.w <- !is.null(w)) {
	    if(n != length(w)) stop("lengths of 'x' and 'w' must match")
	    if(any(w <= 0))    stop("all weights should be positive")
	}
        else w <- rep_len(1, n)
        if(is.unsorted(x)) { # sort along x
            io <- order(x)
            x <- x[io]
            y <- y[io]
            if(has.w)
                w <- w[io]
        }
	Ky <- diff(diff(y)/diff(x))/diff(x,2)
	a0 <- rep(1,n-2)
	out <- .C(cpen, dim = as.integer(c(n,maxit)), t = as.double(x), z = as.double(y), w = as.double(w),
                  a0 = as.double(a0), lambda = as.double(lambda), Ky = as.double(Ky),
                  ## results:
                  L = double(n-1), U = double(n-1),
                  fun = double(n), res = double(n), flag = as.integer(1), tol = as.double(tol), zhat = double(n),
                  iter = as.integer(1), Deriv = double(n))
        ## return
	ret <- list(call = match.call(), x.values = x, y.values = y, lambda = lambda, tol = tol,
                    residuals = out$res, fit.values = out$zhat, alpha = out$a0, deriv = out$Deriv,
                    minvalue = lambda * sum(out$a0 * Ky), lower = out$L, upper = out$U,
                    iter = out$iter, convergence = out$flag, AlphaMVal = out$fun)
	class(ret) <- "cvx.pen.reg"
	ret
}

print.cvx.pen.reg <- function(x, digits = getOption("digits"), ...) {
  cat(sprintf("\"cvx.pen.reg\" object; lambda=%s, tol=%g, Call: \t",
              format(x$lambda, digits=digits), x$tol))
  print(x$call, digits=digits, ...)
  cat(sprintf("Minimum Criterion Value Obtained: %s\n", format(x$minvalue, digits=digits)))
  cat(sprintf("Number of Iterations: %d\n", x$iter))
  cat(sprintf("Convergence Status:   %d\n", x$convergence))
  invisible(x)
}

plot.cvx.pen.reg <- function(x, ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
                             main = sprintf("Convex Regression using\n Penalized Least Squares (lambda=%.4g)", x$lambda),
                             pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...) {
    xx <- x$x.values
    yx <- x$y.values
    fitx <- x$fit.values
    resx <- x$residuals
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0)) -> op; on.exit(par(op))
    plot(xx,yx, xlab = 'x', ylab=ylab, pch=pch, cex=cex, main=main, ...)
    lines(xx, fitx, lwd=lwd, col = col2)
    plot(fitx,resx, xlab = 'Fitted Values', ylab = "Residuals", pch=pch, main = "Fitted vs Residuals", ...)
    abline(h = 0.0, lty = ablty, col = col2)
    plot(yx,fitx, xlab = "Actual Values", ylab = "Fitted Values", pch=pch, main = "Actual vs Fitted", ...)
    abline(a = 0, b = 1, lty = ablty, col = col2)
    qqnorm(resx)
    qqline(resx)
    invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.pen.reg <- function(object, newdata = NULL,...) {
    x <- if(is.null(newdata)) object$x.values else as.vector(newdata)
	## P stores the values of derivatives at newdata.
	## Q stores the function values at newdata.
	## R stores the values of the second derivative at newdata.
	k <- length(x)
	t <- object$x.values
	n <- length(t)
	zhat <- object$fit.values
	deriv <- object$deriv
	L <- object$lower
	U <- object$upper
	fun <- object$AlphaMVal
	out <- .C(predcvxpen, dim = as.integer(c(k,n)), x = as.double(x), t = as.double(t),
                  zhat = as.double(zhat), deriv = as.double(deriv), L = as.double(L), U = as.double(U),
                  fun = as.double(fun), P = double(k), Q = double(k), "2der" = double(k))[
            c("x","Q","P","R")]
    ## return
    cbind(newdata = out$x, fun = out$Q, `1der` = out$P, `2der` = out$R)
}

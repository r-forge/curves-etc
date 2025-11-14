cvx.lse.reg <- function(t, z, w = NULL, ...)
{
    t <- as.vector(t)
    z <- as.vector(z)
    if (!all(is.finite(c(t, z))))
        stop("missing or infinite values in inputs are not allowed")
    if(length(t) != length(z))
        stop("'t' and 'z' must have same length.")
    n <- length(t)
    if(n < 3)
        stop("Number of samples must be at least 3.")
    if(!is.null(w)) {
        if (n != length(w)) stop("lengths of 'x' and 'w' must match")
        if (any(w <= 0))    stop("all weights should be positive")
    }
    if(is.unsorted(t)) { # sort along t; continue to use {t, z, w}, "t = x", "z = y"
        io <- order(t)
        t <- t[io]
        z <- z[io]
        if(!is.null(w))
            w <- w[io]
    }
    ## if(is.null(w)) w <- 1 # may remain scalar!
    n2 <- n - 2L # == nrow(A)
    A <- matrix(0, n2, n) # the tri-diagonal smoother matrix
    for(i in 1:n2) {
        A[i, i:(i+2L)] <- t[c(i+2, i,  i+1)] -
                          t[c(i+1, i+2,i)]
    }
    if(is.null(w)) {
        G <- A
        sw <- w <- 1
    } else { ## use w[]
        sw <- sqrt(w)
        G <- t(t(A)/sw) # ==?  A / rep(sw, each=n-2)
    }
    h <- - A %*% z
    E <- rbind(t(G), t(h))
    f <- c(rep(0,n), 1)
    nnl <- nnls(E,f)
    u <- nnl$x
    r <- as.vector(E %*% u - f) ## FIXME improve '-f' as its mostly -0
    fit <- z + (-r[1:n]/r[n+1])/sw
    deriv <- diff(fit)/diff(t) # length n-1
    deriv <- c(deriv, deriv[n-1L])
    ret <- list(call = match.call(),
                x.values = t, y.values = z, w = w, fit.values = fit, deriv = deriv, iter = 1,
                residuals = z - fit, minvalue = mean(w*(z-fit)^2), convergence = nnl$mode)
    class(ret) <- "cvx.lse.reg"
    ret
}

print.cvx.lse.reg <- function(x, digits = getOption("digits"), ...) {
  cat("Call: \t"); print(x$call)
  cat(sprintf("Minimum Criterion Value Obtained: %s\n", format(x$minvalue, digits=digits)))
  cat(sprintf("Number of Iterations: %d\n", x$iter))
  cat(sprintf("Convergence Status:   %d\n", x$convergence))
  invisible(x)
}

plot.cvx.lse.reg <- function(x, diagnostics = TRUE,
                             ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
                             pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...) {
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  if(diagnostics) {
      plot.window(c(0,7), c(0,7))
      par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0)) -> op; on.exit(par(op))
      plot(xx,yx, xlab = 'x', ylab=ylab, pch=pch, cex=cex,
           main = "Convex Regression using\n Least Squares", ...)
      lines(xx, fitx, lwd=lwd, col = col2)
      plot(fitx,resx, xlab = 'Fitted Values', ylab = "Residuals", pch=pch, main = "Fitted vs Residuals", ...)
      abline(h = 0, lty = ablty, col = col2)
      plot(yx,fitx, xlab = "Actual Values", ylab = "Fitted Values", pch=pch, main = "Actual vs Fitted", ...)
      abline(a = 0, b = 1, lty = ablty, col = col2)
      qqnorm(resx)
      qqline(resx)
  } else {
      plot.window(c(0,7), c(0,7))
      par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.3,.5,0)) -> op; on.exit(par(op))
      plot(xx,yx,xlab = 'x', ylab=ylab, pch=pch, cex=cex, main = "Convex Regression using\n Least Squares", ...)
      lines(xx, fitx, lwd=lwd)
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.lse.reg <- function(object, newdata = NULL, deriv = 0, ...)
{
    stopifnot("deriv must either be 0 or 1!" = (deriv == 0 || deriv == 1))
    zhat <- unname(object$fit.values)
    D <- unname(object$deriv)
    if(is.null(newdata)) {
      warning("No 'newdata' found and so using input 'x' values")
      if(deriv == 0)
          zhat
      else # deriv == 1
          D
    } else {
      newdata <- as.vector(newdata)
      t <- unname(object$x.values)
      dim <- c(length(t), length(newdata), deriv)
             ##  n                r           f
      out <- .C(derivcvxpec, as.integer(dim), as.double(t), as.double(zhat), as.double(D),
                as.double(newdata))
      return(out[[5]])
  }
}

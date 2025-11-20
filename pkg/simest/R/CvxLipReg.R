cvx.lip.reg <- function(t, z, w = NULL, L, ...) {
  t <- as.vector(t)
  z <- as.vector(z)
  if(length(t) != length(z))
      stop("'t' and 'z' must have same length.")
  if (!all(is.finite(c(t, z))))
    stop("missing or infinite values in inputs are not allowed")
  n <- length(t)
  if(n < 3)
      stop("Number of samples must be at least 3.")
  if(!is.null(w)) {
      if (n != length(w)) stop("lengths of 't' and 'w' must match")
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
  A <- matrix(0,nrow = n,ncol = n)
  A[1,] <- c(-1,1,rep(0,n-2))
  for(i in 2:{n-1}){
    A[i,i-1] <- t[i+1] - t[i]
    A[i,i] <- -{t[i+1] - t[i-1]}
    A[i,i+1] <- t[i] - t[i-1]
  }
  A[n,] <- c(rep(0,n-2),1,-1)
  b <- c(-L*(t[2]- t[1]), rep(0,n-2), -L*(t[n] - t[n-1]))
  if(is.null(w)) {
      tG <- t(A)
      sw <- w <- 1
  } else { ## use w[]
      sw <- sqrt(w)
      tG <- t(A)/sw # ==?  A / rep(sw, each=n-2)
  }
  h <- b - A %*% z
  E <- rbind(tG,  t(h))
  f <- c(rep(0,n), 1)
  nnl <- nnls(E, f)
  ##     ====
  tt <- nnl$x
  # print(tt)
  r <- E %*% tt - f ## FIXME improve '-f' as its mostly -0
  fit <- z + (-r[1:n]/r[n+1])/sw
  # print(A%*%fit - b)
  deriv <- diff(fit)/diff(t) # length n-1
  deriv <- c(deriv, deriv[n-1L])
  resi <- z - fit
  ## return
  r <- list(call = match.call(), L = L,
            x.values = t, y.values = z, w = w, fit.values = fit, iter = 1L, deriv = deriv,
            residuals = resi, minvalue = mean(w * resi^2), convergence = nnl$mode)
  class(r) <- "cvx.lip.reg"
  r
}

print.cvx.lip.reg <- function(x, digits = getOption("digits"), ...) {
  cat("Call: \t"); print(x$call)
  cat(sprintf("Minimum Criterion Value Obtained (L=%g): %s\n", x$L, format(x$minvalue, digits=digits)))
  cat(sprintf("Number of Iterations: %d\n", x$iter))
  cat(sprintf("Convergence Status:   %d\n", x$convergence))
  invisible(x)
}

plot.cvx.lip.reg <- function(x, diagnostics = TRUE, ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
                             main = sprintf("Convex Lipschitz Regression\n using Least Squares, L=%g",x$L),
                             pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...) {
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  plot.window(c(0,7), c(0,7))
  par(mfrow = if(diagnostics) c(1,1) else c(2,2),
      mar = c(3,3,3,1), mgp = c(1.3,.5,0)) -> op; on.exit(par(op))
  plot(xx,yx, xlab = 'x', ylab=ylab, pch=pch, cex=cex, main=main, ...)
  lines(xx, fitx, lwd=lwd, col = col2)
  if(diagnostics) { # 3 more plots
      plot(fitx, resx, xlab = 'Fitted Values',ylab = "Residuals",pch=pch,
           main = "Residuals vs Fitted\n(Tukey-Anscombe Plot)", ...)
      abline(h = 0.0, lty = ablty, col = col2)
      plot( yx,  fitx, xlab = "Actual Values",ylab = "Fitted Values",pch=pch,
           main = "Fitted vs Actual", ...)
      abline(a = 0, b = 1, lty = ablty, col = col2)
      qqnorm(resx)
      qqline(resx, col=col2)
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.lip.reg <- function(object, newdata = NULL, deriv = 0, ...) {
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
    	r <- length(xnew <- as.double(newdata))
        t <- unname(object$x.values)
        n <- length(t)
        .C(derivcvxpec, # ../src/derivcvxpec.c
           as.integer(c(n, r, deriv)), as.double(t), as.double(zhat), as.double(D),
           pred = xnew)$pred
    }
}

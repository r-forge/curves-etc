cvx.lip.reg <- function(t, z, w = NULL, L, ...) {
  t <- as.vector(t)
  z <- as.vector(z)
  if(length(t) != length(z))
      stop("'x' and 'y' must have same length.")
  if (!all(is.finite(c(t, z))))
    stop("missing or infinite values in inputs are not allowed")
  n <- length(t)
  if(n < 3)
      stop("Number of samples must be at least 3.")

### MM FIXME: see how to improve in ./CvxLseReg.r -- or better ./CvxLipReg.R.~M++~

  if (is.null(w))
      w <- rep_len(1, n)
  else {
      if(n != length(w)) stop("lengths of 'x' and 'w' must match")
      if(any(w <= 0))    stop("all weights should be positive")
  }
  A <- cbind(t, z, w)
  A <- A[order(A[,1]),]
  x <- as.vector(A[,1])
  y <- as.vector(A[,2])
  w <- as.vector(A[,3])
  A <- matrix(0,nrow = n,ncol = n)
  A[1,] <- c(-1,1,rep(0,n-2))
  for(i in 2:{n-1}){
    A[i,i-1] <- x[i+1] - x[i]
    A[i,i] <- -{x[i+1] - x[i-1]}
    A[i,i+1] <- x[i] - x[i-1]
  }
  A[n,] <- c(rep(0,n-2),1,-1)
  b <- c(-L*(x[2]- x[1]), rep(0,n-2),-L*(x[n] - x[n-1]))
  G <- t(t(A)/sqrt(w))
  h <- b - A%*%y
  E <- rbind(t(G), t(h))
  f <- c(rep(0,n),1)
  tmp <- nnls(E, f)
  tt <- tmp$x
  # print(tt)
  r <- E%*%tt - f
  tt1 <- r[length(r)]
  tt <- -r[-length(r)]/tt1
  fit <- tt/sqrt(w) + y
  # print(A%*%fit - b)
  deriv <- diff(fit)/diff(x)
  deriv <- c(deriv, deriv[length(deriv)])
  ## return
  ret1 <- list(call = match.call(), L = L,
               x.values = x, y.values = y, w = w, fit.values = fit, iter = 1, deriv = deriv,
               residuals = y - fit, minvalue = mean({w*{y-fit}^2}), convergence = tmp$mode)
  class(ret1) <- "cvx.lip.reg"
  ret1
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
      plot(fitx, resx, xlab = 'Fitted Values',ylab = "Residuals",pch=pch,  main = "Fitted vs Residuals", ...)
      abline(h = 0.0, lty = ablty, col = col2)
      plot( yx,  fitx, xlab = "Actual Values",ylab = "Fitted Values",pch=pch,  main = "Actual vs Fitted", ...)
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
    	newdata <- as.vector(newdata)
    	r <- length(newdata)
        t <- unname(object$x.values)
        n <- length(t)
        dim <- c(n, r, deriv)
        out <- .C(derivcvxpec, as.integer(dim), as.double(t), as.double(zhat), as.double(D),
                  as.double(newdata))
        out[[5]]
    }
}

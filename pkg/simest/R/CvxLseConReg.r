cvx.lse.con.reg <- # <<-- to become deprecated
cvx.lse.conreg <- function(t, z, w = NULL, ...)
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
    if(is.null(w)) w <- rep_len(1, n)
    out <- conreg(x=t, y=z, w=w, convex=TRUE, ...)
    ##     ======
    fit <- out$yf
    deriv <- diff(fit)/diff(t)
    deriv <- c(deriv, deriv[length(deriv)])
    ret <- list(call = match.call(),
                 x.values = t, y.values = z, w = w, fit.values = fit, deriv = deriv, iter = 1,
                 residuals = z - fit, minvalue = mean(w*(z-fit)^2), convergence = 1)
    class(ret) <- c("cvx.lse.conreg", "cvx.lse.reg")
                                       #----------> using all the "cvx.lse.reg" methods ---> ./CvxLseReg.r
    ret
}

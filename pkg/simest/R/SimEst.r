sim.est <- function(x, y, w = NULL, beta.init = NULL, nmulti = NULL, L = NULL,
                    lambda = NULL, maxit = 100, bin.tol = 1e-05, beta.tol = 1e-05,
                    method = c("cvx.pen", "cvx.lip", "cvx.lse.con", "cvx.lse", "smooth.pen"),
                    progress = TRUE, force = FALSE)
{
  ## MM {FIXME}: method = "cvs.lse" force = FALSE is deprecated (use  "cvx.lse.con")
  x <- as.matrix(x)
  y <- as.vector(y)
  if (!all(is.finite(c(x,y))))
    stop("missing or infinite values in inputs are not allowed")
  if(nrow(x) != (n <- length(y)))
    stop("Number of rows of 'x' must match the length of 'y'")
  if(n < 3)
    stop("Number of samples must be at least 3.")
  if(!is.null(w)) {
    if(n != length(w)) stop("lengths of 'x' and 'w' must match")
    if(any(w <= 0))    stop("all weights should be positive")
  }
  else w <- rep_len(1,n) # TODO (needed ?) (could keep NULL or use scalar 1)
  method <- match.arg(method)# ==> error if not matching
  ## if (!(method %in% c("cvx.pen", "cvx.lip", "cvx.lse", "smooth.pen","cvx.lse.con")))
  ##   stop("The input for argument 'method' is unrecognized!")
  if (length(lambda) >= 1L && !any(method == c("cvx.pen", "smooth.pen"))) {
    warning("Tuning parameter 'lambda' can only be \nused with 'cvx.pen' or 'smooth.pen'!")
    # method <- "cvx.pen"
  }
  if(length(L) >= 1L && method != "cvx.lip"){
    warning("Tuning parameter 'L' can only used with 'cvx.lip'!")
    # method <- "cvx.lip"
  }
  if (method == "cvx.lse" && !force) # force = FALSE is default  ==> this is _default_ for "cvx.lse"
      method <- "cvx.lse.con"
  else if (method == "cvx.lip" && is.null(L))
    stop("Method 'cvx.lip' needs specification of tuning parameter 'L'")
  else if(method %in% c("cvx.pen","smooth.pen") && is.null(lambda)) {
    lambda <- 0.01*n^{1/5}
    warning(paste0("Required parameter 'lambda' not given; setting lambda = ",lambda))
  }

  beta.path <- function(x, G, tt, xG, GG, tmp) {
    tt4th <- 0.25 * tt * tt * tmp
    ## return
    as.vector((x*(1 + tt4th + tt * xG) - tt * G) / (1 - tt4th))
  }
  smoothFn <-
      switch(method,
             "cvx.pen" = function(t, z, q) {
                 ter <- cvx.pen.reg(t, z, lambda = lambda, w = q)
                 if(min(diff(diff(ter$fit.values)/diff(ter$x.values))) < -1e-01)
                     stop("Function Estimate is not Convex from cvx.pen.reg!!")
                 return(ter)
             },
             "cvx.lip" = function(t, z, q) {
                 ter <- cvx.lip.reg(t, z, w = q, L = L)
                 if(min(diff(diff(ter$fit.values)/diff(ter$x.values))) < -1e-01)
                     stop("Function Estimate is not Convex from cvx.lip.reg!!")
                 return(ter)
             },
             "cvx.lse.con" = function(t, z, q) {
                 ter <- cvx.lse.con.reg(t, z, w = q)
                 if(min(diff(diff(ter$fit.values)/diff(ter$x.values))) < -1e-01)
                     stop("Function Estimate is not Convex from cvx.lse.con.reg!!")
                 return(ter)
             },
             "cvx.lse" = function(t, z, q) {
                 ter <- cvx.lse.reg(t, z, w = q)
                 if(min(diff(diff(ter$fit.values)/diff(ter$x.values))) < -1e-01)
                     stop("Function Estimate is not Convex from cvx.lse.reg!!")
                 return(ter)
             },
             "smooth.pen" = function(t, z, q)
                 smooth.pen.reg(t, z, lambda = lambda, w = q),
             ## otherwise :
             stop("should not happen -- invalid 'method':", method)
             ) # end switch(..)

  d <- ncol(x)
  if(is.null(beta.init)){ # use `nmulti` random beta.init's
    if(is.null(nmulti)){
      nmulti <- round(d*log(d)) + 1
      warning(paste0("'nmulti' unspecified; taking nmulti = ", nmulti))
    }
    beta.init <- matrix(rnorm(nmulti*d), ncol = d)
  } else {
    beta.init <- matrix(beta.init, ncol = d)
    nmulti <- nrow(beta.init)
  }
  beta.init <- beta.init * sign(beta.init[,1]) / sqrt(rowSums(beta.init*beta.init)) # "uniquify": scale to norm 1
  BetaInit <- beta.init
  BetaPath <- matrix(0,nrow = nmulti, ncol = d)
  ObjValPath <- rep_len(1,nmulti)
  FitPath <- vector("list",nmulti)
  xMatPath <- vector("list",nmulti)
  TotIter <- 0L
  ConvCheck <- rep_len(0L, nmulti)
  itervec <-   rep_len(NA_integer_, nmulti)
  for(k in 1:nmulti){
    bfit <- 1e05
    if(progress && k > 1) cat("\rmultistart ", k-1, " of ", nmulti, " done!")
    flag <- 0
    iter <- 0L
    while(iter <= maxit && flag == 0){
      iter <- iter + 1L
      TotIter <- TotIter + 1L
      A <- cbind(x %*% beta.init[k,], y, x)
      tmmp <- fastmerge(A, w = w, tol = bin.tol)
      A <- tmmp$DataMat
      sw <- tmmp$w
      TT <- cbind(sw, A)
      TT <- TT[order(TT[,2]),]
      A <- TT[,-1]; sw <- TT[,1]
      if(min(diff(A[,1])) < 0.9*bin.tol){
        cat("\nMinDiff = ", min(diff(A[,1]))," and bin.tol = ",bin.tol,"\n")
        print(A[,1],20)
        stop("Fastmerge Error!! Datapoints too close!")
      }
      fit <- smoothFn(A[,1], A[,2],sw)

      ## Original {simest} code --------------------------------------------
      G <- -colSums(fit$residuals*fit$deriv*A[, -(1:2), drop=FALSE])
      xG <- sum(beta.init[k,]*G)
      GG <- sum(G*G)
      tmp <- xG*xG - GG
      if(tmp >= 0) { # = 0 not ok as r1 := (...) / (-tmp/2)
        cat("xG*xG - GG = ", tmp,"\n")
        stop("Problem with the Cauchy-Schwarz inequality!")
      }
      bp <- xG - G[1]/beta.init[k,1]
      ## End  Original {simest} code ----------------------------------------

      r1 <- (bp - sqrt(bp*bp - tmp))/(-0.5*tmp)
      r2 <- (bp + sqrt(bp*bp - tmp))/(-0.5*tmp)
      if(is.na(r1) || is.na(r2)) {
          cat("(r1, r2) = ", c(r1, r2),"\n")
          cat(" --- must debug this!!\n") ## later:  stop("............")
          browser()
##############################
      }
      if(bfit - fit$minvalue < beta.tol || iter == maxit){
        if(iter == maxit){
            warning("'maxit' reached for Start no. !!", k, ", iter =", iter)
            ConvCheck[k] <- ConvCheck[k] + 1
        }
        itervec[k] <- iter
        BetaPath[k,] <- beta.init[k,]
        FitPath[[k]] <- fit
        xMatPath[[k]] <- A[, -(1:2), drop=FALSE]
        ObjValPath[k] <- fit$minvalue
        flag <- 1
      } else{
        bfit <- fit$minvalue
      }
      if(r2 - r1 < 2e-05){
        warning("First coordinate stuck at zero.")
        ConvCheck[k] <- ConvCheck[k] + 1
      	itervec[k] <- iter
        BetaPath[k,] <- beta.init[k,]
        FitPath[[k]] <- fit
        xMatPath[[k]] <- A[, -(1:2), drop=FALSE]
        ObjValPath[k] <- fit$minvalue
        flag <- 1
      } else {
        ToOpt <- function(tt){
          A <- cbind(x%*%beta.path(beta.init[k,], G, tt, xG, GG, tmp), y)
          tmmp <- fastmerge(A,w = w,tol=bin.tol)
          A <- tmmp$DataMat
          sw <- tmmp$w
          TT <- cbind(sw, A)
          TT <- TT[order(TT[,2]),]
          A <- TT[,-1]; sw <- TT[,1]
          if(min(diff(sort(A[,1]))) < 0.9*bin.tol){
            cat("MinDiff = ", min(diff(sort(A[,1])))," and bin.tol = ",bin.tol,"\n")
            print(sort(A[,1]),20)
            stop("Fastmerge Error!! Datapoints too close!")
          }
          ffit <- smoothFn(A[,1], A[,2], sw)
          ffit$minvalue
        }
        # opt <- optimize(ToOpt, lower = r1+1e-05, upper = r2-1e-05)$minimum
        opt <- optim(0, ToOpt, method = "L-BFGS-B", lower = r1 + 1e-05, upper = r2 - 1e-05)$par
        beta.init[k,] <- beta.path(beta.init[k,], G, opt, xG, GG, tmp)
      }
    }
  }
  if(progress){
  	cat("\rmultistart ", k, " of ", nmulti, " done!")
  	cat("\n")
  }
  K <- which.min(ObjValPath)
  fit <- FitPath[[K]]
  ## return
  ret <- list(call = match.call(), method = method,
              beta = BetaPath[K,], x.mat = xMatPath[[K]], lambda = lambda, L = L, minvalue = ObjValPath[K],
              nmulti = nmulti, K = K, BetaPath = BetaPath, BetaInit = BetaInit,
              ObjValPath = ObjValPath, regress = fit, iter = TotIter, itervec = itervec,
              x.values = fit$x.values, y.values = fit$y.values, fit.values = fit$fit.values,
              convergence = length(which(ConvCheck == 0)),
              deriv = fit$deriv, residuals = fit$residuals)
  class(ret) <- "sim.est"
  ret
}

print.sim.est <- function(x, digits = getOption("digits"), ...) {
  cat("Call: \t"); print(x$call)
  cat("Estimate of beta is:\n"); print(x$beta, digits=digits, ...)
  cat(sprintf("Number of Starting Vectors: %d\n", x$nmulti))
  cat("Initial vector leading to the minimum:\n")
  print(x$BetaInit[x$K,], digits=digits, ...)
  cat(sprintf("Minimum Criterion Value Obtained: %s\n", format(x$minvalue, digits=digits)))
  cat(sprintf("Number of Iterations: %d\n", x$iter))
  cat(sprintf("Convergence Status:   %d out of %d converged\n", x$convergence, x$nmulti))
  invisible(x)
}

plot.sim.est <- function(x, pch = 20, cex = 1, lwd = 2, col2 = "red", ...) {
    xx <- x$x.values
    yx <- x$y.values
    fitx <- x$fit.values
    resx <- x$residuals
    tt <- switch(x$method,
        "cvx.pen"     = "Convex Link Function Estimate\n using Penalized LS Objective",
        "cvx.lse.con" = "Convex Conreg Link Function Estimate\n using Least Squares Objective",
        "cvx.lse"     = "Convex Link Function Estimate\n using Least Squares Objective",
        "cvx.lip"     = "Convex Link Function Estimate\n using Lipschitz LS Objective",
        "smooth.pen"  = "Unconstrained Link Function Estimate\n using Penalized LS Objective",
        "smooth.pen.gcv" = "Unconstrained Link Function Estimate\n using Penalized LS Objective (GCV)",
        stop(gettextf("invalid 'method' (should not happen, please report!): method=\"%s\"", x$method)))

    plot.window(c(0,7), c(0,7))
    par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(1.8,0.5,0)) -> op; on.exit(par(op))
    plot(xx, yx, xlab = quote(x^T*hat(beta)), ylab = quote(y ~ "and" ~ hat(y) ~ ' values'),
         pch=pch, cex=cex, main = tt, ...)
    lines(xx, fitx, lwd=lwd, col=col2)
    plot(fitx,resx, xlab = 'Fitted Values', ylab = "Residuals", pch=pch, main = "Fitted vs Residuals", ...)
    abline(h = 0, lty = 4, col=col2)
    ## plot(1:x$nmulti, x$ObjValPath,xlab="Start No.",ylab="Objective Value", main = "Minimum Obtained from each Start")
    hist(x$ObjValPath, breaks = x$nmulti, main = "Minimum Obtained from each Start", xlab = "Objective Values")
    qqnorm(resx, main="Normal Q-Q Plot: Residuals", pch=pch)
    qqline(resx)
}

predict.sim.est <- function(object, newdata = NULL, deriv = 0, ...){
    stopifnot("deriv must either be 0 or 1!" = (deriv == 0 || deriv == 1),
              "'newdata' should be numeric!" = 	is.numeric(newdata))
    req <- object$regress
    B <- object$beta
    if(!is.null(newdata)) {
        newdata <- c(newdata)
        newdata <- matrix(newdata, ncol = length(B))
        t <- newdata %*% B
        if(object$method != "cvx.pen")
             predict(req, t, deriv = deriv)
        else predict(req, t)[,deriv+2]
    } else {
        warning("No 'newdata' found and so using input 'x' values")
        if(deriv == 0)
            object$fit.values
        else # (deriv == 1)
            object$deriv
    }
}

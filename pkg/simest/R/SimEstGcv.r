simestgcv <- function(x, y, w = NULL, beta.init = NULL, nmulti = NULL,
                      lambda = NULL, maxit = 100, bin.tol = 1e-06, beta.tol = 1e-05,
                      agcv.iter = 100, progress = TRUE)

{
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- length(y)
  if(n <= 2)
    stop("Number of samples must be greater than 2.")
  if (!all(is.finite(cbind(x, y))))
    stop("missing or infinite values in inputs are not allowed")
  if(nrow(x) != length(y))
    stop("Number of rows of 'x' must match the length of 'y'")
  if(is.null(lambda))
    stop("Required input 'lambda' missing.")
  if(length(lambda) == 2)
    lambda <- c(min(lambda), max(lambda))
  beta.path <- function(x, G, tt, xG, GG, tmp){
    tt4th <- 0.25 * tt * tt * tmp
    ## return
    as.vector((x*(1 + tt4th + tt * xG) - tt * G) / (1 - tt4th))
  }
  smoothFn <- function(t, z, q, lam) smooth.pen.reg(t, z, lambda = lam, w = q)
  d <- ncol(x)
  if(is.null(w))
    w <- rep_len(1,n)
  else {
    if(length(w) != n)
      stop("'w' and 'y' should be of same length!")
    if(any(w <= 0))
      stop("'w' should contain positive elements!")
  }
  if(is.null(beta.init)){
    if(is.null(nmulti)){
      nmulti <- round(d*log(d)) + 1
      warning(paste("'nmulti' not given. \nTaking nmulti = ",nmulti,sep = ""))
    }
    beta.init <- matrix(rnorm(nmulti*d),ncol = d)
  } else{
    beta.init <- matrix(beta.init,ncol = d)
    nmulti <- nrow(beta.init)
  }
  beta.init <- beta.init * sign(beta.init[,1]) / sqrt(rowSums(beta.init^2))

  GcvOpt <- function(lam, BetaInit) {
    ObjValPath <- rep_len(1,nmulti)

    chkmergeA <- function(A, msg) { # using global {bin.tol, lam}
        if((md <- min(diff(A1 <- A[,1]))) < 0.9*bin.tol) { # bail out (1)
            cat("MinDiff and bin.tol are:","\n")
            print(c(MinD = md, bin.tol=bin.tol), digits=16)
            print(sort(A1), digits=16)
            stop("fastmerge() error in GcvOpt(lam=%g) -- Datapoints too close! - in ", lam, msg, call.=FALSE)
        }
    }
    for(k in 1:nmulti) {
      bfit <- 1e5 # current "best fit" obj.fn. value
      flag <- 0
      iter <- 0
      while(iter <= maxit && flag == 0) {
        iter <- iter + 1
        A <- cbind(x %*% BetaInit[k,], y, x)
        tmmp <- fastmerge(A, w = w, tol = bin.tol)
        A <- tmmp$DataMat
        sw <- tmmp$w
        TT <- cbind(sw, A)
        TT <- TT[order(TT[,2]),]
        A <- TT[,-1]; sw <- TT[,1]
        chkmergeA(A, msg = sprintf("1: k=%d, iter=%d", k,iter))
        fit <- smoothFn(A[,1], A[,2],sw,lam)
        ##     ~~~~~~~~
        G <- -colSums(fit$residuals*fit$deriv*A[, -(1:2), drop=FALSE])
        xG <- sum(BetaInit[k,]*G)
        GG <- sum(G*G)
        tmp <- xG*xG - GG
        if(tmp > 0) {
            stop("xG*xG - GG = ", tmp," > 0: something went wrong!")
        }
        bp <- xG - G[1]/BetaInit[k,1]
        r1 <- (bp - sqrt(bp*bp - tmp))/(-0.5*tmp)
        r2 <- (bp + sqrt(bp*bp - tmp))/(-0.5*tmp)
        if(is.na(r1) || is.na(r2)) cat("(r1, r2) = ", c(r1, r2),"\n")
        if(bfit - fit$minvalue < beta.tol || iter == maxit) {
          if(iter == maxit){
            warning("'maxit' reached for start no.", k, " out of ", nmulti)
          }
          ObjValPath[k] <- fit$minvalue
          flag <- 1
        } else{
          bfit <- fit$minvalue
        }
        if(r2 - r1 < 2e-05){
          warning("First coordinate stuck at zero.")
          ObjValPath[k] <- fit$minvalue
          flag <- 1
        } else {
          ToOpt <- function(tt) {
            A <- cbind(x %*% beta.path(BetaInit[k,], G, tt, xG, GG, tmp), y)
            tmmp <- fastmerge(A, w = w, tol=bin.tol)
            A <- tmmp$DataMat
            sw <- tmmp$w
            TT <- cbind(sw, A)
            TT <- TT[order(TT[,2]),]
            A <- TT[,-1]; sw <- TT[,1]
            chkmergeA(A, msg = sprintf("2: inside ToOpt(tt=%g), k=%d, iter=%d", tt, k,iter))
            ffit <- smoothFn(A[,1],A[,2],sw,lam)
            ffit$minvalue
          }
          # opt <- optimize(ToOpt, lower = r1+1e-05, upper = r2-1e-05)$minimum
          opt <- optim(0, ToOpt, method = "L-BFGS-B", lower = r1 + 1e-05, upper = r2 - 1e-05)$par
          BetaInit[k,] <- beta.path(BetaInit[k,], G, opt, xG, GG, tmp)
        }
      } # end while(iter <= maxit ...)
      ##
    } # end for(in 1:nmulti)
    ##  --- ---
    K <- which.min(ObjValPath)
    A <- cbind(x %*% BetaInit[K,], y, x)
    tmmp <- fastmerge(A, w = w, tol = bin.tol)
    A <- tmmp$DataMat
    sw <- tmmp$w
    TT <- cbind(sw, A)
    TT <- TT[order(TT[,2]),]
    A <- TT[,-1]; sw <- TT[,1]
    chkmergeA(A, msg = sprintf("finally before last smooth.pen.reg(): iter=%d, K=%d", iter, K))
    ## return the computed  agcv.score :
    smooth.pen.reg(A[,1], A[,2], w=sw, lambda=lam, agcv=TRUE, agcv.iter=agcv.iter)$agcv.score
  } # end{ GcvOpt }

  ## Find lambda := arg min_{\lambda} GcvOpt(\lambda, \beta_0, ...)
  ##
  ## BestGcv <- optim(0.01*n^{1/5}, GcvOpt, method = "L-BFGS-B", lower = lambda[1],
  		# upper = lambda[2], BetaInit = beta.init)
  BestGcv <- optimize(GcvOpt, lower = lambda[1], upper = lambda[2], BetaInit = beta.init, maximum  = FALSE)
  lam <- BestGcv$minimum
  ret <- sim.est(x = x,y = y, method = "smooth.pen", lambda = lam, beta.init = beta.init, progress = progress)
  ##     =======
  ret$ method <- "smooth.pen.gcv"
  ret$ GCVoptim <- BestGcv
  ret$ call <- match.call()
  ret
}

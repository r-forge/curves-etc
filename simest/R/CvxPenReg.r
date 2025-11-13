cvx.pen.reg <- function(t, z, lambda, w = NULL, max.iter = 500, alpha.tol= 1e-03, ...) UseMethod("cvx.pen.reg")

cvx.pen.reg.default <- function(t, z, lambda, w = NULL, max.iter = 500, alpha.tol= 1e-03, ...){
  p <- lambda
  tol <- alpha.tol
  if(length(t) != length(z))
    stop("'x' and 'y' must have same length.")  
  if (!all(is.finite(c(t, z)))) 
    stop("missing or infinite values in inputs are not allowed")
  n <- length(t)
  if(n <= 2)
    stop("Number of samples must be greater than 2.")
  w <- if (is.null(w)) 
        rep_len(1, n)/n
    else {
        if (n != length(w)) 
            stop("lengths of 'x' and 'w' must match")
        if (any(w <= 0)) 
            stop("all weights should be positive")
        (w * sum(w > 0))/n
    }
  q <- w
  if (!is.finite(p) || p < 0)
    stop("'lambda' must be non-negative and finite")  
  A <- cbind(t, z, w)
  A <- A[order(A[,1]),]
  t <- A[,1]
  z <- A[,2]
  w <- A[,3]
  n <- length(t)
  N <- n-2
  K <- matrix(0,nrow = N, ncol = N+2)
  a.i <- diff(diff(t))/4
  for(i in 1:N){
    K[i,i] <- 1/{{t[i+1]-t[i]}*{t[i+2]-t[i]}}
    K[i,i+1] <- -1/{{t[i+2]-t[i+1]}*{t[i+1]-t[i]}}
    K[i,i+2] <- 1/{{t[i+2]-t[i]}*{t[i+2]-t[i+1]}}
  }
  a.init <- a.i
  Mfun <- function(a){
    mf <- .C("Mfun", as.integer(N), as.double(a), as.double(t), double(2*N + 2), 
    	double(N*N), PACKAGE = "simest")
    Vmat <- mf[[4]]
    V <- matrix(Vmat, ncol = 2, nrow = N+1, byrow = TRUE)
    Mmat <- mf[[5]]
    M <- matrix(Mmat, ncol = N, nrow = N, byrow = TRUE)
    return(list(M = M, V = V))
  }
  norm = function(a){
    return(sqrt(sum(a^2)))
  }
  M11 <- p*(K%*%diag(1/q))%*%t(K); M22 <- K%*%z
  ba <- a.init
  MM <- Mfun(ba)$M
  bf <- norm({{MM + M11}%*%ba} - M22)
  run <- function(a0, sup.iter = max.iter){
    loop <- rep(0,sup.iter/50)
    for(i in 1:{sup.iter}){
      MF <- Mfun(a0)
      Bf0 <- MF$M + M11
      a1 <- solve.pentadiag(Bf0, M22)
      check <- norm(Bf0%*%a0 - M22)
      if(check < bf){
        ba <- a0
        bf <- check
      }
      if(i%%50 == 0){
        if(is.element(check,loop)){
          #print(paste('Code Looped at stop criterion value',check))
          #warning('Fixed point iteration started looping.\n And so returning the best value obtained.')
          flag <- 2
          temp <- as.vector(p*diag(1/q)%*%t(K)%*%ba)
          temp1 <- t(ba)%*%M22
          zhat <- z - temp
          V <- MF$V
          DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(ba), 
          	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
          D <- DD[[length(DD)]]
          return(list(alpha = ba, minvalue = c(p*temp1), x.values = t, y.values = z, 
            fit.values = zhat, residuals = temp, convergence = flag, 
            Kz = M22, Vmat = MF$V, deriv = D, iter = i))
        }
        else loop[i/50] <- check
      }
      if(check <= tol){
        temp <- as.vector(p*diag(1/q)%*%t(K)%*%a1)
        temp1 <- t(a1)%*%M22
        zhat <- z - temp
        flag <- 0
        V <- Mfun(a1)$V
        DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(a1),
        	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
        D <- DD[[length(DD)]]
        return(list(alpha = a1, minvalue = c(p*temp1), x.values = t, y.values = z, 
          fit.values = zhat, residuals = temp, convergence = flag,
          Kz = M22, Vmat = MF$V, deriv = D, iter = i))
      }
      if(i == sup.iter){
        #warning('sup.iter in convex function estimation reached.')
        flag <- 1
        temp <- as.vector(p*diag(1/q)%*%t(K)%*%ba)
        temp1 <- t(ba)%*%M22
        zhat <- z - temp
        V <- MF$V
        DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(ba),
        	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
        D <- DD[[length(DD)]]
        return(list(alpha = a1, minvalue = c(p*temp1), x.values = t, y.values = z, 
          fit.values = zhat, residuals = temp, convergence = flag, 
          Kz = M22, Vmat = MF$V, deriv = D, iter = i))
      }
      else a0 <- a1
      i <- i + 1
    }
  }
  ret <- run(a.init)
  ret$call <- match.call()
  class(ret) <- "cvx.pen.reg"
  return(ret)
}

print.cvx.pen.reg <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Minimum Criterion Value Obtained:\n")
  print(x$minvalue)
  cat("Number of Iterations:\n")
  print(x$iter)
  cat("Convergence Status:\n")
  print(x$convergence)  
  #cat("Input 'x', 'y', the fit and the derivative:\n")
  #write.table(format(rbind(c('x','y','yhat','deriv'),
  #  cbind(object$x.values, object$y.values, object$fit.values, object$deriv)), 
  #  justify="right"), row.names=F, col.names = F, quote=F)
}

plot.cvx.pen.reg <- function(x, ...){
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  diagnostics = TRUE
  if(diagnostics){
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
    plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Penalized Least Squares")
    lines(xx, fitx, lwd = 2)
    plot(fitx,resx,xlab = 'Fitted Values',ylab = "Residuals",pch = "*", type = 'p', main = "Fitted vs Residuals")
    abline(h = 0.0, lty = 4)
    plot(yx,fitx,xlab = "Actual Values",ylab = "Fitted Values",pch = "*", type = 'p', main = "Actual vs Fitted")
    abline(a = 0, b = 1, lty = 4)
    qqnorm(resx)
    qqline(resx)
  } else{
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
    plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Penalized Least Squares")
    lines(xx, fitx, lwd = 2)    
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.pen.reg <- function(object,newdata = NULL,...){
  t <- object$x.values
  zhat <- object$fit.values
  a <- object$alpha
  D <- object$deriv
  V <- object$Vmat
  n <- length(t)
  if(is.null(newdata)){
    warning("No 'newdata' found and so using input 'x' values")
    return(zhat)
  } else{
    newdata <- as.vector(newdata)
    r <- length(newdata)
    foo <- function(kk){
      if(kk == t[1])
        return(zhat[1])
      if(kk < t[1])
        return(zhat[1] + D[1]*{kk - t[1]})
      if(kk > t[n])
        return(zhat[n] + D[n]*{kk - t[n]})
      if(kk > t[1] & kk <= t[2]){
        return({a[1]*(a[1] > 0)*{kk - t[1]}^3/{6*{t[2] - t[1]}*{t[3] - t[1]}}} 
          + D[1]*{kk - t[1]} + zhat[1])
      }
      if(kk > t[n-1] & kk <= t[n]){
        ret <- {a[n-2]*(a[n-2] > 0)*{kk - t[n-1]}*{t[n] - t[n-1]}/{2*{t[n] - t[n-2]}}} + {a[n-2]*(a[n-2] > 0)*{{t[n] - kk}^3 - {t[n] - t[n-1]}^3}/{6*{t[n] - t[n-1]}*{t[n] - t[n-2]}}} + {D[n-1]*{kk - t[n-1]}} + zhat[n-1]
        return(ret)
      }
      for(i in 3:{n-1}){
        if(kk > t[i-1] & kk <= t[i]){
          if(kk <= V[i-1,2]){
            trm1 <- -a[i-2] * {kk - t[i]}^3 / {6*{t[i] - t[i-1]} * {t[i] - t[i-2]}}
            trm2 <- -a[i-2]*{t[i-1] - t[i]}^2/{6*{t[i] - t[i-2]}}
            trm3 <- a[i-1]*{kk - t[i-1]}^3/{6*{t[i]-t[i-1]}*{t[i+1] - t[i-1]}}
            trm4 <- a[i-2]*{kk - t[i-1]}*{t[i] - V[i-1,1]}^2/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
            trm5 <- -a[i-1]*{kk - t[i-1]}*{V[i-1,1] - t[i-1]}^2/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}
            trm6 <- D[i-1]*{kk - t[i-1]}
            trm7 <- zhat[i-1]
            ret <- trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7
            return(ret)
          }else{
            kk1 <- V[i-1,2]
            trm1 <- -a[i-2]*{kk1 - t[i]}^3/{6*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
            trm2 <- -a[i-2]*{t[i-1] - t[i]}^2/{6*{t[i] - t[i-2]}}
            trm3 <- a[i-1]*{kk1 - t[i-1]}^3/{6*{t[i]-t[i-1]}*{t[i+1] - t[i-1]}}
            trm4 <- a[i-2]*{kk1 - t[i-1]}*{t[i] - V[i-1,1]}^2/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
            trm5 <- -a[i-1]*{kk1 - t[i-1]}*{V[i-1,1] - t[i-1]}^2/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}
            trm6 <- D[i-1]*{kk1 - t[i-1]}
            trm7 <- zhat[i-1]
            ret <- trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7
            ret <- ret
                 + D[i-1]*{kk - kk1}
                 + {{a[i-2]*{{t[i] - V[i-1,1]}^2 - {t[i] - V[i-1,2]}^2}/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}} + {a[i-1]*{{V[i-1,2] - t[i-1]}^2 - {V[i-1,1] - t[i-1]}^2}/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}}}*{kk - V[i-1,2]}
            return(ret)
          }
        }
      }
    }
    return(as.vector(unlist(sapply(newdata[seq_len(r)],foo,simplify = TRUE,USE.NAMES = FALSE))))
  }
}
# funest <- function(t,zhat,a,D,V,n){
#   foo <- function(kk){
#     if(kk == t[1]){
#       return(zhat[1])
#     }
#     if(kk < t[1]){
#       ret <- zhat[1] + D[1]*{kk - t[1]}
#       return(ret)
#     }
#     if(kk > t[n]){
#       ret <- zhat[n] + D[n]*{kk - t[n]}
#       return(ret)
#     }
#     if(kk > t[1] & kk <= t[2]){
#       ret <- {a[1]*(a[1] > 0)*{kk - t[1]}^3/{6*{t[2] - t[1]}*{t[3] - t[1]}}} + D[1]*{kk - t[1]} + zhat[1]
#       return(ret)
#     }
#     if(kk > t[n-1] & kk <= t[n]){
#       ret <- {a[n-2]*(a[n-2] > 0)*{kk - t[n-1]}*{t[n] - t[n-1]}/{2*{t[n] - t[n-2]}}} + {a[n-2]*(a[n-2] > 0)*{{t[n] - kk}^3 - {t[n] - t[n-1]}^3}/{6*{t[n] - t[n-1]}*{t[n] - t[n-2]}}} + {D[n-1]*{kk - t[n-1]}} + zhat[n-1]
#       return(ret)
#     }
#     for(i in 3:{n-1}){
#       if(kk > t[i-1] & kk <= t[i]){
#         if(kk <= V[i-1,2]){
#           trm1 <- -a[i-2] * {kk - t[i]}^3 / {6*{t[i] - t[i-1]} * {t[i] - t[i-2]}}
#           trm2 <- -a[i-2]*{t[i-1] - t[i]}^2/{6*{t[i] - t[i-2]}}
#           trm3 <- a[i-1]*{kk - t[i-1]}^3/{6*{t[i]-t[i-1]}*{t[i+1] - t[i-1]}}
#           trm4 <- a[i-2]*{kk - t[i-1]}*{t[i] - V[i-1,1]}^2/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
#           trm5 <- -a[i-1]*{kk - t[i-1]}*{V[i-1,1] - t[i-1]}^2/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}
#           trm6 <- D[i-1]*{kk - t[i-1]}
#           trm7 <- zhat[i-1]
#           ret <- trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7
#           return(ret)
#         }else{
#           kk1 <- V[i-1,2]
#           trm1 <- -a[i-2]*{kk1 - t[i]}^3/{6*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
#           trm2 <- -a[i-2]*{t[i-1] - t[i]}^2/{6*{t[i] - t[i-2]}}
#           trm3 <- a[i-1]*{kk1 - t[i-1]}^3/{6*{t[i]-t[i-1]}*{t[i+1] - t[i-1]}}
#           trm4 <- a[i-2]*{kk1 - t[i-1]}*{t[i] - V[i-1,1]}^2/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}
#           trm5 <- -a[i-1]*{kk1 - t[i-1]}*{V[i-1,1] - t[i-1]}^2/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}
#           trm6 <- D[i-1]*{kk1 - t[i-1]}
#           trm7 <- zhat[i-1]
#           ret <- trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7
#           ret <- ret
#                + D[i-1]*{kk - kk1}
#                + {{a[i-2]*{{t[i] - V[i-1,1]}^2 - {t[i] - V[i-1,2]}^2}/{2*{t[i] - t[i-1]}*{t[i] - t[i-2]}}} + {a[i-1]*{{V[i-1,2] - t[i-1]}^2 - {V[i-1,1] - t[i-1]}^2}/{2*{t[i] - t[i-1]}*{t[i+1] - t[i-1]}}}}*{kk - V[i-1,2]}
#           return(ret)
#         }
#       }
#     }
#   }
#   return(foo)
# }
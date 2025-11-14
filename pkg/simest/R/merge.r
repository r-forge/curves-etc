fastmerge <- function(DataMat, w = NULL, tol = 1e-04){
  DataMat <- as.matrix(DataMat)
  p <- ncol(DataMat)
  x <- DataMat[,1]
  n <- length(x)
  if(is.null(w)) w <- rep(1,n)
  sigT <- function(x) floor(x) + tol*floor((x - floor(x))/tol)
  xx <- sigT(x)
  ## clearly this has been based on  smooth.spline()'s way ...
  nd <- !duplicated(xx)
  ux <- sort(x[nd])
  uxx <- sort(xx[nd])
  nx <- length(ux)
  if (nx == n) {
    ox <- TRUE
    wDmat <- cbind(w, DataMat, 0)
  } else {
    ox <- match(xx, uxx)
## FIXME use vapply()
    tapply1 <- function(X, INDEX, FUN = NULL, ..., simplify = TRUE){
      sapply(X = unname(split(X, INDEX)), FUN = FUN, ..., simplify = simplify, USE.NAMES = FALSE)
    }
    meanV <- function(i, D, q)
        c(sum(q[i]), colMeans(D[i,,drop = FALSE]),
          if((ni <- length(i)) == 1L) 0 else var(D[i,2])*(ni-1L))
    wDmat <- matrix(unlist(use.names = FALSE,
                           tapply1(seq_len(n), ox, meanV, D = DataMat, q = w)),
                    ncol = p+2, byrow = TRUE)
  }
  w <- wDmat[, 1L]
  DataMat <- wDmat[, -1L]
  DataMat[,1L] <- sigT(DataMat[,1L]) # 'x'
  ## return
  list(DataMat = DataMat[,-(p+1L)], w = w, AddVar = DataMat[,p+1L])
}

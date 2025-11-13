solve.pentadiag <- function(a, b,...){
	a <- as.matrix(a)
	if(nrow(a) != ncol(a))
		stop("'a' is not a square matrix!")
	b <- as.vector(b)
	if(length(b) != ncol(a))
		stop("'a' and 'b' should be of same length!")
	n <- length(b)
	aaa <- diag(a)
	bb <- c(diag(a[-n,-1]), 0)
	c <- c(diag(a[-c(n, n-1), -c(1,2)]), 0, 0)
	d <- c(0, diag(a[-1,-n]))
	e <- c(0, 0, diag(a[-c(1,2), -c(n,n-1)]))
	u <- v <- w <- p <- q <- r <- s <- t <- aa <- dd <- rep_len(0,n+2)
	z <- rep_len(0,n)
	BB <- .C("penta", as.integer(n), as.double(aaa), as.double(bb), as.double(c), as.double(d),
	as.double(e), as.double(b), as.double(u), as.double(v), as.double(w), as.double(p), as.double(q),
	as.double(r), as.double(s), as.double(t), as.double(aa), as.double(dd), as.double(z), 
	PACKAGE = "simest")
	z <- as.vector(BB[[length(BB)]])
	return(z)
}
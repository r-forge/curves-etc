solve_pentadiag <- function(a, b) {
	a <- as.matrix(a)
	if(!is.double(a)) storage.mode(a) <- "double" # --> all other are double
	if((n <- nrow(a)) != ncol(a)) stop("'a' is not a square matrix!")
        if(n < 3) stop("Dimension n must be at least 3")
	b <- as.vector(b)
	if(length(b) != n) stop("'b' must have 'length(b) == nrow(a)'")
	D <- diag(a)
        in1 <- seq_len(n - 1L)
        in2 <- in1[-(n - 1L)] # 1:(n-2)
	C <- c(a[cbind(in1, i2n <- 2:n)], 0)
	F <- c(a[cbind(in2, i3n <- 3:n)], 0, 0)
	A <- c(a[cbind(i2n, in1)], 0)
	E <- c(a[cbind(i3n, in2)], 0, 0)
        ## return    //  C code in ../src/penta.c
        .C(penta, as.integer(n), E, A, D, C, F,
           as.double(b),  z = double(n))$z
}

solve.pentadiag <- function(a, b) {
    .Deprecated("solve_pentadiag")
    solve_pentadiag(a, b)
}

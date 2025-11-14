\name{fastmerge}
\title{Pre-binning of Data Points}
\alias{fastmerge}
\description{
  Numerical tolerance problems in non-parametric regression makes it
  necessary for pre-binning of data points. This procedure is implicitly
  performed by most of the regression functions in \R. This function
  implements this procedure with a given tolerance level.
}
\usage{
fastmerge(DataMat, w = NULL, tol = 1e-4)
}
\arguments{
  \item{DataMat}{a numeric matrix/vector with rows as data points.}
  \item{w}{an optional numeric vector of the same length as \eqn{x}; defaults to all elements 1.}
  \item{tol}{a numeric value providing the tolerance for identifying
    duplicates with respect to the first column \code{DataMat[,1]}.}
}

\details{
  If two values in the first column of \code{DataMat} are separated by a
  value less than \code{tol} then the corresponding rows are merged.
}
\value{
  A list including the elements
  \item{DataMat}{a numeric matrix/vector with rows sorted and possibly
    merged with respect to the first column.}
  \item{w}{obtained weights corresponding to the merged points.}
}
\seealso{
  The function \code{\link{smooth.spline}} also uses such pre-binning.
}
\author{Arun Kumar Kuchibhotla;
  also the authors of \code{\link{smooth.spline}}.
}
\examples{ ## relevant example % found in ../tests/fastmerge-ex.R
n <- 47
set.seed(2657) # <- found after searching
x <- sort(signif(runif(n, -1,1), 5))
y <- sinpi(3*x) * exp(-x) + rnorm(n)/10
str(fmL <- fastmerge(cbind(x,y))) # only 44 (out of 47) "unique" x[]
d.fm <- data.frame(fmL)
d2 <- data.frame(fastmerge(cbind(x,y), tol = 25e-4)) # larger tol ==> only 42 "unique"
table(w <- d2$w) # 3x w=2  and  1 w=3
stopifnot(nrow(d.fm) == 44, nrow(d2) == 42,
          identical(    w[w > 1], c(2, 2, 2, 3)),
          identical(which(w > 1), c(5L, 26L, 28L, 39L)),
          all.equal(1000 * fmL$AddVar[fmL$w != 1],
                    c(2.28919, 23.918, 17.5813), tolerance = 3e-6))

plot(y ~ x, type = "b")
lines(d.fm[,1], d.fm[,2], col = adjustcolor(2, 1/2), lwd=3)
lines(d2  [,1], d2  [,2], col = adjustcolor(4, 1/2), lwd=2)
abline(v = d.fm[d.fm$w > 1, 1],  col = 2, lwd=3, lty=2)
abline(v = (xw <- d2[w > 1, 1]), col = 4, lwd=2, lty=3)
axis(3, at= xw, labels=paste("w=",w[w > 1]), col = 4, col.axis = 4)
}
\concept{Pre-binning}

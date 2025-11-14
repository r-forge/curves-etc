\name{cvx.lse.con.reg}
\title{Convex Least Squares Regression}
\alias{cvx.lse.conreg}
\alias{cvx.lse.con.reg}% to be deprecated
\description{
 Provides an estimate of the non-parametric regression
  function with a shape constraint of convexity and no smoothness
  constraint.  Note that convexity by itself provides some implicit
  smoothness.
}
\usage{
cvx.lse.conreg(t, z, w = NULL, ...)
cvx.lse.con.reg(t, z, w = NULL, ...) # will be deprecated
}
\arguments{
  \item{t, z}{numeric vectors (of the same lengths) with the values of the
    predictor and response variable.}
  \item{w}{an optional numeric vector of the same length as \code{t};
    defaults to all equal weights.}% all equal \eqn{1/n}.
  \item{...}{additional arguments, passed to \code{\link[cobs]{conreg}()}.}
}
\details{
  This function does the same thing as \code{cvx.lse.reg} except that here
  we use the \code{\link[cobs]{conreg}} function from \CRANpkg{cobs}
  package which is faster than \code{cvx.lse.reg}.

  The \code{\link{plot}}, \code{\link{predict}}, \code{\link{print}} functions of \code{\link{cvx.lse.reg}}
  also apply for \code{cvx.lse.con.reg}.
}
\value{
  An object of class \code{cvx.lse.reg}, basically a list including the elements
  \item{x.values}{sorted \code{t} values provided as input.}
  \item{y.values}{corresponding \code{z} values in input.}
  \item{fit.values}{corresponding fit values of same length as that of \code{x.values}.}
  \item{deriv}{corresponding values of the derivative of same length as that of \code{x.values}.}
  \item{iter}{number of steps taken to complete the iterations.}
  \item{residuals}{residuals obtained from the fit.}
  \item{minvalue}{minimum value of the objective function attained.}
  \item{convergence}{a numeric indicating the convergence of the code. Always set to 1.}
}
\source{
  Lawson, C. L and Hanson, R. J. (1995)
  \emph{Solving Least Squares Problems.} SIAM.
}
\references{
  Chen, D. and Plemmons, R. J. (2009)
  \emph{Non-negativity Constraints in Numerical Analysis.}
  Symposium on the Birth of Numerical Analysis.

  Liao, X. and Meyer, M. C. (2014).
  \code{coneproj}: An R package for the primal or dual cone projections with routines for constrained regression.
  \emph{Journal of Statistical Software} \bold{61}(12), 1--22.
}
\author{Arun Kumar Kuchibhotla}% no longer working: arunku@wharton.upenn.edu
\examples{
str(cvx.lse.conreg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
tmp <- cvx.lse.con.reg(x, y)
print(tmp)
plot(tmp)
predict(tmp, newdata = rnorm(10,0,0.1))
}
\keyword{regression}
\keyword{nonlinear}
\keyword{optimize}
\concept{Convex Least Squares}
\concept{Cone Projection}

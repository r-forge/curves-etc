\name{cvx.lse.reg}
\title{Convex Least Squares Regression}
\alias{cvx.lse.reg}
\alias{plot.cvx.lse.reg}
\alias{predict.cvx.lse.reg}
\alias{print.cvx.lse.reg}
\description{
  This function provides an estimate of the non-parametric regression
  function with a shape constraint of convexity and no smoothness
  constraint. Note that convexity by itself provides some implicit
  smoothness.
}
\usage{
cvx.lse.reg(t, z, w = NULL,...)
\method{plot}{cvx.lse.reg}(x, diagnostics = TRUE,
     ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
     pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...)
\method{print}{cvx.lse.reg}(x, digits = getOption("digits"), ...)
\method{predict}{cvx.lse.reg}(object, newdata = NULL, deriv = 0, ...)
}

\arguments{
  \item{t}{a numeric vector giving the values of the predictor variable.}
  \item{z}{a numeric vector giving the values of the response variable.}
  \item{w}{an optional numeric vector of the same length as t; Defaults to all elements \eqn{1/n}.}
  \item{...}{additional arguments.}

  \item{diagnostics}{for the \code{\link{plot}()} method; if true, as by
    default, produce diagnostics, notably residual plots additionally.}
  \item{ylab, pch, cex, lwd, col2, ablty}{further optional argument to the
    \code{plot()} method; the last two for the color and line type of
    \emph{some} plot components.}
  \item{digits}{the number of significant digits, for numbers in the
    \code{print()} method.}

  \item{x, object}{an object of class \code{"cvx.lse.reg"}.}
  \item{newdata}{a matrix of new data points in the predict function.}
  \item{deriv}{a numeric either 0 or 1 representing which derivative to evaluate.}
}
\details{
  The function minimizes \deqn{\sum_{i=1}^n w_i(z_i - \theta_i)^2,}
  subject to
  \deqn{\frac{\theta_2 - \theta_1}{t_2 - t_1}\le\cdots\le\frac{\theta_n - \theta_{n-1}}{t_n - t_{n-1}},}
  %%
  for sorted \eqn{t} values (and \eqn{z} permuted accordingly such that \eqn{(t_i, z_i)} stay pairs.

  This function previously used the \code{coneA()} function from the
  \CRANpkg{coneproj} package to perform the constrained minimization of
  least squares. Currently, the code makes use
  of the \code{\link[nnls]{nnls}()} function from package \CRANpkg{nnls}
  for the same purpose.

  The \code{plot} method provides a scatterplot along with the fitted
  curve; it also includes some diagnostic plots for residuals.
  The \code{\link{predict}()} method allows computation of the first
  derivative.
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
  \item{convergence}{a numeric indicating the convergence of the code.}
}
\source{
  Lawson, C. L and Hanson, R. J. (1995)
  \emph{Solving Least Squares Problems.} SIAM.
}
\references{
  Chen, D. and Plemmons, R. J. (2009)
  \emph{Non-negativity Constraints in Numerical Analysis.}
  Symposium on the Birth of Numerical Analysis.

  Liao, X. and Meyer, M. C. (2014)
  \code{coneproj}: An R package for the primal or dual cone projections with routines for constrained regression.
  \emph{Journal of Statistical Software} \bold{61}(12), 1--22.
}
\author{Arun Kumar Kuchibhotla}% no longer arunku@wharton.upenn.edu}
\examples{
args(cvx.lse.reg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
cvxL <- cvx.lse.reg(x, y)
print(cvxL)
plot(cvxL)
predict(cvxL, newdata = rnorm(10,0,0.1))
}
\concept{Convex Least Squares}
\concept{Cone Projection}

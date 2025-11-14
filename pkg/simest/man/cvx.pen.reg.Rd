\name{cvx.pen.reg}
\title{Penalized Smooth Convex Regression}
\alias{cvx.pen.reg}
\alias{plot.cvx.pen.reg}
\alias{predict.cvx.pen.reg}
\alias{print.cvx.pen.reg}
\description{
  This function provides an estimate of the non-parametric regression
  function with a shape constraint of convexity and smoothness constraint
  provided through square integral of second derivative.
}
\usage{
cvx.pen.reg(x, y, lambda, w = NULL, tol = 1e-5, maxit = 1000)

\method{plot}{cvx.pen.reg}(x, ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
     main = sprintf("Convex Regression using\n Penalized Least Squares (lambda=\%.4g)",
                    x$lambda),
     pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...)
\method{print}{cvx.pen.reg}(x, digits = getOption("digits"), ...)
\method{predict}{cvx.pen.reg}(object, newdata = NULL,...)
}

\arguments{
  \item{x}{a numeric vector giving the values of the predictor variable.
    For \code{plot} and \code{print} functions, \code{x} is an object of
    class \code{"cvx.pen.reg"}.}
  \item{y}{a numeric vector giving the values of the response variable.}
  \item{lambda}{a nonnegative number, the penalty value aka \sQuote{smoothing parameter}.}
  \item{w}{an optional numeric vector of the same length as x; defaults to all 1.}
  \item{tol}{non-negative number, the tolerance level for convergence.}
  \item{maxit}{an integer giving the maxmimum number of steps taken by the algorithm;
    defaults to 1000.}
  %---
  \item{ylab, main, pch, cex, lwd, col2, ablty}{further optional argument to the
    \code{plot()} method; the last two for the color and line type of
    \emph{some} plot components.}
  \item{...}{any additional arguments, passed, e.g., to lower level plotting functions.}
  \item{digits}{the number of significant digits, for numbers in the
    \code{print()} method.}
  \item{object}{an object of class \code{"cvx.pen.reg"}; for the \code{\link{predict}()} method.}
  \item{newdata}{a vector of new data points to be used in \code{predict()}.}
}

\details{
  The function minimizes
  \deqn{\sum_{i=1}^n w_i(y_i - f(x_i))^2 + \lambda\int\{f''(x)\}^2dx}
  subject to convexity constraint on \eqn{f}.

  \code{plot} function provides the scatterplot along with fitted curve; it
  also includes some diagnostic plots for residuals. Predict function
  returns a matrix containing the inputted newdata along with the function
  values, derivatives and second derivatives.
}
\value{
  An object of class \code{cvx.pen.reg}, basically a list including the elements
  \item{x.values}{sorted \code{x} values provided as input.}
  \item{y.values}{corresponding \code{y} values in input.}
  \item{fit.values}{corresponding fit values of same length as that of \code{x.values}.}
  \item{deriv}{corresponding values of the derivative of same length as that of \code{x.values}.}
  \item{iter}{number of steps taken to complete the iterations.}
  \item{residuals}{residuals obtained from the fit.}
  \item{minvalue}{minimum value of the objective function attained.}
  \item{convergence}{a numeric indicating the convergence of the code.}
  \item{alpha}{a numeric vector of length 2 less than \code{x}.  This
    represents the coefficients of the B-splines in the second derivative
    of the estimator.}
  \item{AlphaMVal}{a numeric vector needed for predict function.}
  \item{lower}{a numeric vector needed for predict function.}
  \item{upper}{a numeric vector needed for predict function.}
}
\source{
  Elfving, T. and Andersson, L. (1988).
  An Algorithm for Computing Constrained Smoothing Spline Functions.
  \emph{Numer. Math.} \bold{52}(5), 583--595.

  Dontchev, A. L., Qi, H. and Qi, L. (2003).
  Quadratic Convergence of Newton's Method for Convex Interpolation and
  Smoothing.
  \emph{Constructive Approximation} \bold{19}(1),123--143.
}
\author{Arun Kumar Kuchibhotla, arunku@wharton.upenn.edu,
  Rohit Kumar Patra, rohit@stat.columbia.edu.}
\examples{
args(cvx.pen.reg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
tmp <- cvx.pen.reg(x, y, lambda = 0.01)
print(tmp)
plot(tmp)
predict(tmp, newdata = rnorm(10,0,0.1))
}
\keyword{Penalized Least Squares}

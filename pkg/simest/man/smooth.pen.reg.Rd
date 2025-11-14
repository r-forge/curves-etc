\name{smooth.pen.reg}
\title{Penalized Smooth/Smoothing Spline Regression}
\alias{smooth.pen.reg}
\alias{plot.smooth.pen.reg}
\alias{predict.smooth.pen.reg}
\alias{print.smooth.pen.reg}
%% \docType{data}
\description{
  Estimate the non-parameteric regression function using smoothing splines.
}
\usage{
smooth.pen.reg(x, y, lambda, w = NULL, agcv = FALSE, agcv.iter = 100, ...)

\method{plot}{smooth.pen.reg}(x, diagnostics = TRUE,
     ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
     pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...)
\method{print}{smooth.pen.reg}(x, digits = getOption("digits"), ...)
\method{predict}{smooth.pen.reg}(object, newdata = NULL, deriv = 0, ...)
}

\arguments{
  \item{x}{a numeric vector giving the values of the predictor variable.
    For plot and print methods, \code{x} is an object of class \code{smooth.pen.reg}.}
  \item{y}{a numeric vector giving the values of the response variable.}
  \item{lambda}{a numeric value giving the penalty value.}
  \item{w}{an optional numeric vector of the same length as \code{x}; Defaults to all 1.}
  \item{agcv}{a logical denoting if an estimate of generalized cross-validation is needed.}
  \item{agcv.iter}{a numeric denoting the number of random vectors used to
    estimate the GCV.  See \sQuote{Details}.}

  \item{...}{additional arguments, passed further, e.g., to
    \code{\link{matplot}} or \code{\link{plot.default}} for the \code{plot()} method.}

  \item{diagnostics}{for the \code{\link{plot}()} method; if true, as by
    default, produce diagnostics, notably residual plots additionally.}
  \item{ylab, pch, cex, lwd, col2, ablty}{further optional argument to the
    \code{plot()} method; the last two for the color and line type of
    \emph{some} plot components.}
  \item{digits}{the number of significant digits, for numbers in the
    \code{print()} method.}
  \item{object}{the result of \code{smooth.pen.reg()}, of class \code{smooth.pen.reg}.}
  \item{newdata}{a matrix of new data points for the \code{\link{predict}} method.}
  \item{deriv}{either 0 or 1, the order of the derivative to evaluate.}
}

\details{
  The function minimizes
  \deqn{\sum_{i=1}^n w_i(y_i - f(x_i))^2 + \lambda\int\{f''(x)\}^2\; dx}
  without any constraint on \eqn{f}.

  This function implements in \R the algorithm noted in Green and
  Silverman(1994).  The function \code{\link{smooth.spline}} in \R is not
  suitable for single index model estimation as it chooses \eqn{\lambda} using GCV by
  default.

  \code{plot} function provides the scatterplot along with fitted
  curve; it also includes some diagnostic plots for residuals. Predict
  function now allows computation of the first derivative. Calculation of
  generalized cross-validation requires the computation of diagonal
  elements of the hat matrix involved which is cumbersone and is
  computationally expensive (and also is unstable).

  \code{smooth.Pspline()} in the \CRANpkg{pspline} package provides the GCV
  criterion value which matches the usual GCV when all the weights are
  equal to 1 but is not clear what it is for weights unequal. We use an
  estimate of GCV (formula of which is given in Green and Silverman (1994))
  proposed by Girard which is very stable and computationally cheap. For
  more details about this randomized GCV, see Girard (1989).
}
\value{
  An object of class \code{smooth.pen.reg}, basically a list including the elements
  \item{x.values}{sorted \code{x} values provided as input.}
  \item{y.values}{corresponding \code{y} values in input.}
  \item{fit.values}{corresponding fit values of same length as that of \code{x.values}.}
  \item{deriv}{corresponding values of the derivative of same length as that of \code{x.values}.}
  \item{iter}{Always set to 1.}
  \item{residuals}{residuals obtained from the fit.}
  \item{minvalue}{minimum value of the objective function attained.}
  \item{convergence}{Always set to 0.}
  \item{agcv.score}{Asymptotic GCV approximation. Proposed in Silverman (1982) as a computationally fast approximation to GCV.}
  \item{splinefun}{An object of class \code{smooth.spline} needed for predict.}
}
\references{
  Green, P. J. and Silverman, B. W. (1994)
  \emph{Non-parametric Regression and Generalized Linear Models: A Roughness Penalty Approach.}
  Chapman and Hall.

  Girard, D. A. (1989)
  A Fast Monte-Carlo Cross-Validation Procedure for Large Least Squares
  Problems with Noisy Data.
  \emph{Numerische Mathematik} \bold{56}, 1--23.
}
\author{Arun Kumar Kuchibhotla}
\examples{
args(smooth.pen.reg)
x <- runif(50,-1,1)
y <- x^2 + 0.3 * rnorm(50)
smP <- smooth.pen.reg(x, y, lambda = 0.01, agcv = TRUE)
smP # -> print() method
plot(smP)
predict(smP, newdata = rnorm(10, 0,0.1))
}
\keyword{Smoothing Splines}
\keyword{Penalized Least Squares}

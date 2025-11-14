\name{cvx.lip.reg}
\title{Convex Least Squares Regression with Lipschitz Constraint}
\alias{cvx.lip.reg}
\alias{plot.cvx.lip.reg}
\alias{predict.cvx.lip.reg}
\alias{print.cvx.lip.reg}
\description{
  This function provides an estimate of the non-parametric regression
  function with a shape constraint of convexity and a smoothness constraint
  via a Lipschitz bound.
}
\usage{
cvx.lip.reg(t, z, w = NULL, L, ...)

\method{plot}{cvx.lip.reg}(x, diagnostics = TRUE, ylab = quote(y ~ "and" ~ hat(y) ~ " values"),
     main = sprintf("Convex Lipschitz Regression\n using Least Squares, L=\%g", x$L),
     pch = "*", cex = 1, lwd = 2, col2 = "red", ablty = 4, ...)
\method{print}{cvx.lip.reg}(x, digits = getOption("digits"), ...)
\method{predict}{cvx.lip.reg}(object, newdata = NULL, deriv = 0, ...)
}
\arguments{
  \item{t}{a numeric vector giving the values of the predictor variable.}
  \item{z}{a numeric vector giving the values of the response variable.}
  \item{w}{an optional numeric vector of the same length as \code{t}; Defaults to all elements \eqn{1/n}.}
  \item{L}{a numeric value providing the Lipschitz bound on the function.}
  \item{...}{additional arguments.}

  \item{diagnostics}{for the \code{\link{plot}()} method; if true, as by
    default, produce diagnostic, notably residual plots additionally.}
  \item{main, ylab, pch, cex, lwd, col2, ablty}{further optional argument
    to the \code{plot()} method; the last two for the color and line type
    of \emph{some} plot components.}
  \item{digits}{the number of significant digits, for numbers in the
    \code{print()} method.}
  \item{x, object}{an object of class \code{"cvx.lip.reg"}.}
  \item{newdata}{a matrix of new data points in the predict function.}
  \item{deriv}{a numeric either 0 or 1 representing which derivative to evaluate.}
}

\details{
  The function minimizes  \deqn{\sum_{i=1}^n w_i(z_i - \theta_i)^2}
  subject to
  \deqn{-L \le \frac{\theta_2-\theta_1}{t_2-t_1} \le \cdots \le \frac{\theta_n-\theta_{n-1}}{t_n-t_{n-1}} \le L,}{
        -L <=     (\theta_2-\theta_1) / (t_2-t_1) <=   ...    <=   (\theta_n-\theta_{n-1}) / (t_n-t_{n-1}) <=   L,}
  for sorted \eqn{t} values (and \eqn{z} permuted accordingly such that \eqn{(t_i, z_i)} stay pairs.

  This function uses the \code{nnls} function from the \CRANpkg{nnls} package to
  perform the constrained minimization of least squares.  The
  \code{\link{plot}} method provides the scatterplot along with fitted
  curve; when \code{diagnostics = TRUE}, it also includes some diagnostic
  plots for residuals.

  The \code{\link{predict}()} method allows calculating the first
  derivative as well.
}
 \value{
  An object of class \code{"cvx.lip.reg"}, basically a list with elements
  \item{x.values}{sorted \code{t} values provided as input.}
  \item{y.values}{corresponding \code{z} values in input.}
  \item{fit.values}{corresponding fit values of same length as that of \code{x.values}.}
  \item{deriv}{corresponding values of the derivative of same length as that of \code{x.values}.}
  \item{residuals}{residuals obtained from the fit.}
  \item{minvalue}{minimum value of the objective function attained.}
  \item{iter}{always set to 1, here.}
  \item{convergence}{a numeric indicating the convergence status of the code.}
}
\source{
  Lawson, C. L and Hanson, R. J. (1995).
  Solving Least Squares Problems. SIAM.
}
\seealso{
  Function \code{\link[nnls]{nnls}} from CRAN package \CRANpkg{nnls}.
}
\references{
  Chen, D. and Plemmons, R. J. (2009).
  Non-negativity Constraints in Numerical Analysis.
  Symposium on the Birth of Numerical Analysis.
}
\author{Arun Kumar Kuchibhotla}% no longer working: arunku@wharton.upenn.edu
\examples{
args(cvx.lip.reg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
tmp <- cvx.lip.reg(x, y, L = 10)
print(tmp)
plot(tmp)
predict(tmp, newdata = rnorm(10,0,0.1))
}
\keyword{regression}
\keyword{nonlinear}
\keyword{optimize}
\concept{Convex Least Squares}
\concept{Least Distance Programming}
\concept{Non-negative Least Squares}

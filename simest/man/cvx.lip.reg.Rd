\name{cvx.lip.reg}
\alias{cvx.lip.reg}
\alias{cvx.lip.reg.default}
\alias{plot.cvx.lip.reg}
\alias{predict.cvx.lip.reg}
\alias{print.cvx.lip.reg}
\docType{data}
\title{
Convex Least Squares Regression with Lipschitz Constraint
}
\description{
This function provides an estimate of the non-parametric regression function with a shape constraint of convexity and a smoothness constraint via a Lipschitz bound.
}
\usage{
cvx.lip.reg(t, z, w = NULL, L,...)
\method{cvx.lip.reg}{default}(t, z, w = NULL, L, ...)
\method{plot}{cvx.lip.reg}(x,...)
\method{print}{cvx.lip.reg}(x,...)
\method{predict}{cvx.lip.reg}(object, newdata = NULL,...)
}

\arguments{
	\item{t}{a numeric vector giving the values of the predictor variable.}
	\item{z}{a numeric vector giving the values of the response variable.}
	\item{w}{an optional numeric vector of the same length as x; Defaults to all elements \eqn{1/n}.}
	\item{L}{a numeric value providing the Lipschitz bound on the function.}
	\item{...}{additional arguments.}
	\item{x}{an object of class `cvx.lip.reg'.}
	\item{object}{An object of class `cvx.lip.reg'.}
	\item{newdata}{a matrix of new data points in the predict function.}
}

\details{
The function minimizes 
\deqn{\sum_{i=1}^n w_i(y_i - \theta_i)^2} 
subject to
\deqn{-L\le\frac{\theta_2 - \theta_1}{x_2 - x_1}\le\cdots\le\frac{\theta_n - \theta_{n-1}}{x_n - x_{n-1}}\le L}
for sorted \eqn{x} values and \eqn{y} reorganized such that \eqn{y_i} corresponds to the new sorted \eqn{x_i}. This function uses the \code{nnls} function from the \code{nnls} package to perform the constrained minimization of least squares. \code{plot} function provides the scatterplot along with fitted curve; it also includes some diagnostic plots for residuals.
}
\value{
	An object of class `cvx.lip.reg', basically a list including the elements
	\item{x.values}{sorted `x' values provided as input.}
	\item{y.values}{corresponding `y' values in input.}
	\item{fit.values}{corresponding fit values of same length as that of `x.values'.}
	\item{deriv}{corresponding values of the derivative of same length as that of `x.values'.}
	\item{residuals}{residuals obtained from the fit.}
	\item{minvalue}{minimum value of the objective function attained.}
	\item{iter}{Always set to 1.}
	\item{convergence}{a numeric indicating the convergence of the code.}
}
\source{
Lawson, C. L and Hanson, R. J. (1995). Solving Least Squares Problems. SIAM.
}
\seealso{
	See also the function \code{\link{nnls}}.
}
\references{
Chen, D. and Plemmons, R. J. (2009). Non-negativity Constraints in Numerical Analysis. Symposium on the Birth of Numerical Analysis. 
}
\author{Arun Kumar Kuchibhotla, arunku@wharton.upenn.edu.}
\examples{
args(cvx.lip.reg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
tmp <- cvx.lip.reg(x, y, L = 10)
print(tmp)
plot(tmp)
predict(tmp, newdata = rnorm(10,0,0.1))
}
\keyword{Convex Least Squares}
\keyword{Least Distance Programming}
\keyword{Non-negative Least Squares}
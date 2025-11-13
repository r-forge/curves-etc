\name{cvx.lse.reg}
\alias{cvx.lse.reg}
\alias{cvx.lse.reg.default}
\alias{plot.cvx.lse.reg}
\alias{predict.cvx.lse.reg}
\alias{print.cvx.lse.reg}
\docType{data}
\title{
Convex Least Squares Regression.
}
\description{
This function provides an estimate of the non-parametric regression function with a shape constraint of convexity and no smoothness constraint. Note that convexity by itself provides some implicit smoothness.
}
\usage{
cvx.lse.reg(t, z, w = NULL,...)
\method{cvx.lse.reg}{default}(t, z, w = NULL, ...)
\method{plot}{cvx.lse.reg}(x,...)
\method{print}{cvx.lse.reg}(x,...)
\method{predict}{cvx.lse.reg}(object, newdata = NULL,...)
}

\arguments{
	\item{t}{a numeric vector giving the values of the predictor variable.}
	\item{z}{a numeric vector giving the values of the response variable.}
	\item{w}{an optional numeric vector of the same length as t; Defaults to all elements \eqn{1/n}.}
	\item{...}{additional arguments.}	
	\item{x}{An object of class `cvx.lse.reg'. This is for plot and print function.}
	\item{object}{An object of class `cvx.lse.reg'.}
	\item{newdata}{a matrix of new data points in the predict function.}
}
\details{
The function minimizes 
\deqn{\sum_{i=1}^n w_i(y_i - \theta_i)^2} 
subject to
\deqn{\frac{\theta_2 - \theta_1}{x_2 - x_1}\le\cdots\le\frac{\theta_n - \theta_{n-1}}{x_n - x_{n-1}}}
for sorted \eqn{x} values and \eqn{y} reorganized such that \eqn{y_i} corresponds to the new sorted \eqn{x_i}. This function uses the \code{coneA} function from the \code{coneproj} package to perform the constrained minimization of least squares.
}
\value{
	An object of class `cvx.lse.reg', basically a list including the elements
	\item{x.values}{sorted `x' values provided as input.}
	\item{y.values}{corresponding `y' values in input.}
	\item{fit.values}{corresponding fit values of same length as that of `x.values'.}
	\item{deriv}{corresponding values of the derivative of same length as that of `x.values'.}
	\item{iter}{number of steps taken to complete the iterations.}
	\item{residuals}{residuals obtained from the fit.}
	\item{minvalue}{minimum value of the objective function attained.}
	\item{convergence}{a numeric indicating the convergence of the code.}
}
\source{
Liao, X. and Meyer, M. C. (2014). coneproj: An R package for the primal or dual cone projections with routines for constrained regression. Journal of Statistical Software 61(12), 1 -- 22.
}
\author{Arun Kumar Kuchibhotla, arunku@wharton.upenn.edu}
\examples{
args(cvx.lse.reg)
x <- runif(50,-1,1)
y <- x^2 + rnorm(50,0,0.3)
tmp <- cvx.lse.reg(x, y)
print(tmp)
plot(tmp)
predict(tmp, newdata = rnorm(10,0,0.1))
}
\keyword{Convex Least Squares}
\keyword{Cone Projection}
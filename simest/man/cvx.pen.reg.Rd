\name{cvx.pen.reg}
\alias{cvx.pen.reg}
\alias{cvx.pen.reg.default}
\alias{plot.cvx.pen.reg}
\alias{predict.cvx.pen.reg}
\alias{print.cvx.pen.reg}
\docType{data}
\title{
Penalized Smooth Convex Regression.
}
\description{
This function provides an estimate of the non-parametric regression function with a shape constraint of convexity and smoothness constraint provided through square integral of second derivative.
}
\usage{
cvx.pen.reg(t, z, lambda, w = NULL, 
	max.iter = 500, alpha.tol = 1e-03,...)
\method{cvx.pen.reg}{default}(t, z, lambda, w = NULL, 
	max.iter = 500, alpha.tol = 1e-03, ...)
\method{plot}{cvx.pen.reg}(x,...)
\method{print}{cvx.pen.reg}(x,...)
\method{predict}{cvx.pen.reg}(object, newdata = NULL,...)
}

\arguments{
	\item{t}{a numeric vector giving the values of the predictor variable.}
	\item{z}{a numeric vector giving the values of the response variable.}
	\item{lambda}{a numeric value giving the penalty value.}
	\item{w}{an optional numeric vector of the same length as x; Defaults to all 1.}
	\item{max.iter}{an integer giving the maxmimum number of steps taken by the algorithm; defaults to 500.}
	\item{alpha.tol}{a numeric providing the tolerance level for convergence.}
	\item{...}{additional arguments.}
	\item{x}{An object of class `cvx.pen.reg'. This is for plot and print function.}
	\item{object}{An object of class `cvx.pen.reg'.}
	\item{newdata}{a matrix of new data points in the predict function.}
}

\details{
The function minimizes 
\deqn{\sum_{i=1}^n w_i(y_i - f(x_i))^2 + \lambda\int\{f''(x)\}^2dx} 
subject to convexity constraint on \eqn{f}.
}
\value{
	An object of class `cvx.pen.reg', basically a list including the elements
	\item{x.values}{sorted `x' values provided as input.}
	\item{y.values}{corresponding `y' values in input.}
	\item{fit.values}{corresponding fit values of same length as that of `x.values'.}
	\item{deriv}{corresponding values of the derivative of same length as that of `x.values'.}
	\item{iter}{number of steps taken to complete the iterations.}
	\item{residuals}{residuals obtained from the fit.}
	\item{minvalue}{minimum value of the objective function attained.}
	\item{convergence}{a numeric indicating the convergence of the code.}
	\item{alpha}{a numeric vector of length 2 less than `x'. A vector needed for predict.}
	\item{Kz}{a numeric vector of length 2 less than `x' giving second order divided difference of \eqn{y} with respect to \eqn{x}. Needed for predict.}
	\item{Vmat}{a numeric matrix needed for predict function.}
}
\source{
Elfving, T. and Andersson, L. (1988). An Algorithm for Computing Constrained Smoothing Spline Functions. Numer. Math., 52(5):583--595.
}
\author{Arun Kumar Kuchibhotla, arunku@wharton.upenn.edu.}
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
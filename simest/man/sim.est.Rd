\name{sim.est}
\alias{sim.est}
\alias{sim.est.default}
\alias{print.sim.est}
\alias{plot.sim.est}
\alias{predict.sim.est}
\docType{data}
\title{
Single Index Model Estimation: Objective Function Approach.
}
\description{
This function provides an estimate of the non-parametric function and the index vector by minimizing an objective function specified by the method argument.
}
\usage{
sim.est(x, y, method = c("cvx.pen","cvx.lip","cvx.lse","smooth.pen"),
		lambda = NULL, beta.init = NULL, w = NULL, nmulti = 5, bin.tol = 1e-04,
		beta.tol = 1e-04, beta.iter = 100, msg = TRUE,...)
\method{sim.est}{default}(x, y, method = c("cvx.pen","cvx.lip","cvx.lse","smooth.pen"),
		lambda = NULL, beta.init = NULL, w = NULL, nmulti = 5, bin.tol = 1e-04,
		beta.tol = 1e-04, beta.iter = 100, msg = TRUE,...)
\method{plot}{sim.est}(x,...)
\method{print}{sim.est}(x,...)
\method{predict}{sim.est}(object, newdata = NULL,...)
}

\arguments{
	\item{x}{a numeric matrix giving the values of the predictor variables or covariates. For functions plot and print, `x' is an object of class `sim.est'.}
	\item{y}{a numeric vector giving the values of the response variable.}
	\item{method}{a string indicating which method to use for regression.}
	\item{lambda}{a numeric value giving the penalty value.}
	\item{w}{an optional numeric vector of the same length as \eqn{x}; Defaults to all 1.}
	\item{beta.init}{An numeric vector giving the initial value for the index vector.}
	\item{nmulti}{An integer giving the number of multiple starts to be used for iterative algorithm. If beta.init is provided then the nmulti is set to 1.}
	\item{bin.tol}{A tolerance level upto which the x values used in regression are recognized as distinct values.}
	\item{beta.tol}{A tolerance level for stopping iterative algorithm for the index vector.}
	\item{beta.iter}{Maximum number of iterations for the index vector.}
	\item{msg}{If msg is true then some messages will be printed while the code executes.}
	\item{...}{Any additional arguments.}
	\item{object}{An object of class `sim.est'.}
	\item{newdata}{a matrix of new data points in the predict function.}
}

\details{
The function minimizes 
\deqn{\sum_{i=1}^n w_i(y_i - f(x_i^{\top}\beta))^2 + \lambda\int\{f''(x)\}^2dx}
with constraints on \eqn{f} dictated by method. The penalty parameter \eqn{\lambda} is not choosen by any criteria. It has to be specified for using method = `cvx.pen.reg' or `smooth.pen.reg' and \eqn{L} has to be specified for using method \eqn{=} `cvx.lip.reg'. 
}
\value{
	An object of class `smooth.pen.reg', basically a list including the elements
	\item{beta}{A numeric vector storing the estimate of the index vector.}
	\item{nmulti}{Number of multistarts used.}
	\item{x.mat}{the input `x' matrix with possibly aggregated rows.}
	\item{BetaPath}{a list containing the paths taken by each initial index vector for nmulti times.}
	\item{ObjValPath}{a matrix with nmulti rows storing the path of objective function value for multiple starts.}
	\item{convergence}{a vector storing convergence status for the index parameter and the function estimation.}
	\item{itervec}{a vector of length nmulti storing the number of iterations taken by each of the multiple starts.}
	\item{iter}{a numeric giving the total number of iterations taken.}
	\item{method}{method given as input.}
	\item{regress}{An output of the regression function used needed for predict.}
	\item{x.values}{sorted `x.betahat' values obtained by the algorithm.}
	\item{y.values}{corresponding `y' values in input.}
	\item{fit.values}{corresponding fit values of same length as that of `x.betahat'.}
	\item{deriv}{corresponding values of the derivative of same length as that of `x.betahat'.}
	\item{residuals}{residuals obtained from the fit.}
	\item{minvalue}{minimum value of the objective function attained.}
}
\source{
Kuchibhotla, A. K., Patra, R. K. and Sen, B. (2015+). On Single Index Models with Convex Link.
}
\author{Arun Kumar Kuchibhotla, arunku@wharton.upenn.edu}
\examples{
args(sim.est)
x <- matrix(runif(50*2,-1,1),ncol = 2)
b0 <- rep_len(1,2)/sqrt(2)
y <- (x\%*\%b0)^2 + rnorm(50,0,0.3)
#tmp1 <- sim.est(x, y, lambda = 0.01, method = "cvx.pen")
tmp2 <- sim.est(x, y, lambda = 10, method = "cvx.lip")
#tmp3 <- sim.est(x, y, lambda = 0.01, method = "smooth.pen")
tmp4 <- sim.est(x, y, method = "cvx.lse")
#print(tmp1)
print(tmp2)
#print(tmp3)
print(tmp4)
#plot(tmp1)
plot(tmp2)
#plot(tmp3)
plot(tmp4)
#predict(tmp1, newdata = c(0,0,0))
predict(tmp2, newdata = c(0,0,0))
#predict(tmp3, newdata = c(0,0,0))
predict(tmp4, newdata = c(0,0,0))
}
\keyword{Single Index Model}
\keyword{Penalized Least Squares}
\title{lpridge}{
Local polynomial regression fitting with ridging
}
\usage{
lpridge(x, y, bandwidth, deriv=0, n.out=200, x.out=NULL,
	order=NULL, ridge=NULL, weight="epa", mnew=100,
	var=F)
}
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as x.}
  \item{bandwidth}{bandwidth for nonparametric estimation.  Either a
    number or a vector of the same length as x.out. }
  \item{deriv}{order of derivative of the regression function to be
    estimated; default is 0.}
  \item{n.out}{number of output design points, where the function has to
    be estimated.  The default value is n.out=200.}
  \item{x.out}{vector of output design points, where the function has to
    be estimated.  The default value is an equidistant grid of n.out
    points from min(x) to max(x).}
  \item{order}{order of the polynomial used for local polynomials. The
    default value is order=deriv+1.}
  \item{ridge}{ridging parameter. The default value performs a slight
    ridging (see DETAILS below). ridge=0 leads to the local polynomial
    estimator without ridging.}
  \item{weight}{kernel weight function.  The default value is
    weight="epa" for Epanechnikov weights.  Other weights are "bi" for
    biweights (square of "epa") and "tri" for triweights (cube of
    "epa"). If weight is a vector, it is interpreted as vector of
    coefficients of the polynomial weight function.  Thus, weight="epa"
    is equivalent to weight=c(1,0,-1).}
  \item{mnew}{force to restart the algorithm after mnew updating
    steps. The default value is mnew=100. For mnew=1 you get a
    numerically "super-stable" algorithm (see ref. SBE\\&G below).}
  \item{var}{logical flag: if TRUE, the variance of the estimator
    proportional to the residual variance is computed (see DETAILS
    below).}
}
\value{
  a list including used parameters and estimator.
  \item{x}{vector of ordered design points.}
  \item{y}{vector of observations ordered according to x.}
  \item{bandwidth}{vector of bandwidths actually used for nonparametric
    estimation.}
  \item{deriv}{order of derivative of the regression function estimated.}
  \item{x.out}{vector of ordered output design points.}
  \item{order}{order of the polynomial used for local polynomials.}
  \item{ridge}{ridging parameter used.}
  \item{weight}{vector of coefficients of the kernel weight function.}
  \item{mnew}{force to restart the algorithm after mnew updating steps. }
  \item{var}{logical flag: whether the variance of the estimator was computed.}
  \item{est}{estimator of the derivative of order deriv of the
    regression function.}
  \item{est.var}{estimator of the variance of est (proportional to
    residual variance).}
}
\description{
  Fast and stable algorithm for nonparametric estimation of regression
  functions and their derivatives via local polynomials and local
  polynomial ridge regression with polynomial weight functions.
}
\details{described in the reference SBE\\&G below. Several parameters
  described in SBE\\&G are fixed either in the fortran routine or in the
  Splus-function. There, you find comments how to change them.

  In S\\&G a bad finite sample behavior of local polynomials for random
  design was found, and ridging of the estimator was proposed. In this
  algorithm, we use a ridging matrix corresponding to the smoothness
  assumption \emph{``The squared difference of the derivative of order
    deriv of the regression function at the point of estimation and the
    weighted mean of design points is bounded by the residual variance
    divided by the ridge parameter.''}

  Thus, without any smoothness assumption, ridge=0 would be appropriate,
  and for a nearly constant derivative of order deriv, a ridge parameter
  going to infinity behaves well. For equidistant design, ridging
  influences the estimator only at the boundary. Asymptotically, the
  influence of any non-increasing ridge parameter vanishes.

  So far, our experience with the choice of a ridging parameter is
  limited. Therefore we have chosen a default value, which performs a
  slight modification of the local polynomial estimator (with
  denotations h=bandwidth, d=deriv, and where
  n_0=length(x)*mean(bandwidth)/diff(range(x)) is a mean number of
  observations in a smoothing interval):

  ridge = 5*sqrt(n_0)*h^(2*d)/((2*d+3)*(2*d+5))

  For var=T, the variance of the estimator proportional to the residual
  variance is computed, i.e. the exact finite sample variance of the
  regression estimator is var(est) = est.var*sigma^2.
}
\references{The same as for \code{\link{lpepa}}.
}
\examples{
myfit <- lpridge(x,y,b,ridge=0)		# local polynomials
plot(x,y)
lines(myfit$x.out,myfit$est,col=2)
myridge <- lpridge(x,y,b)		# local pol. ridge
lines(myridge$x.out,myridge$est,col=3)
}
\keyword{smoothing}
% Converted by Sd2Rd version 0.2-a1.
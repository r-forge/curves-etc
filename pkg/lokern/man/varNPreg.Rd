\name{varNPreg}
\alias{varNPreg}
\alias{varest}
\title{Nonparametric Variance Estimator}
\description{
  Estimates the error variance \eqn{\sigma^2} nonparametrically in the model
  \deqn{Y_i = m(x_i) + E_i,} where
  \eqn{E_i \sim (0,\sigma^2)}{E_i ~ (0,sigma^2)}, i.i.d.

  Computes leave-one-out residuals (local linear approximation followed by
  reweighting) and their variance.
}
\usage{
varNPreg(x, y)
}
\arguments{
  \item{x}{abscissae values, ordered increasingly.}
  \item{y}{observations at \code{y[i]} at \code{x[i]}.}
}
\note{
  This is an \R interface to the \code{resest} Fortran subroutine, used
  in \code{\link{lokerns}} and \code{\link{glkerns}}, see the latter's help
  page for references and context.

  Earlier version of the \pkg{lokern} package accidentally contained
  \code{varest()} which has been an identical copy of \code{varNPreg()}.
}
\value{
  A list with components
  \item{res}{numeric; residuals at \code{x[]} of length \code{n}.}
  \item{snr}{explained variance of the true curve, i.e., an \eqn{R^2},
    defined as \eqn{1 - \hat{\sigma^2}/ \hat{\sigma_0^2}}{1 - sigma^2/sigm0^2}, where
    \eqn{\hat{\sigma^2} = }{sigma^2 =}\code{sigma2}, and
    \eqn{\hat{\sigma_0^2} := var(Y) = E[Y^2] - (E[Y])^2}{sigm0^2 = var(Y) = E[Y^2] - (E[Y])^2},
    see the example below.} % -> ../src/auxkerns.f
  \item{sigma2}{estimation of residual variance, \eqn{\hat{\sigma^2}}{\sigma^2}.}
}
\author{Martin Maechler}
\seealso{\code{\link{lokerns}}, \code{\link{glkerns}}.}
\examples{
n <- 100
x <- sort(runif(n))
y <- sin(pi*x) + rnorm(n)/10
str(ve <- varNPreg(x,y))

plot(x, y)
## "fitted" = y - residuals:
lines(x, y - ve$res, col=adjustcolor(2, 1/2), lwd=3)
segments(x,y,x,y-ve$res, col=3:4, lty=2:3, lwd=1:2)

## sigma2 :=  1/n  sum_i res_i^2  :
          with(ve,         c(sigma2, sum(res^2)/n))
stopifnot(with(ve, all.equal(sigma2, sum(res^2)/n)))

## show how 'snr' is computed, given 'sigma2'  { in ../src/auxkerns.f }
dx2 <- diff(x, 2) # (x[i+1] - x[i-1]) i= 2..{n-1}
dx.n <- c(x[2]-x[1], dx2, x[n]-x[n-1])
SY  <- sum(dx.n * y)
SY2 <- sum(dx.n * y^2)
rx <- 2*(x[n]-x[1]) # 'dn'
(sigm2.0 <- SY2/rx - (SY/rx)^2)
(R2 <- 1 - ve$sigma2 / sigm2.0)
stopifnot(all.equal(ve$snr, R2))
}
\keyword{nonparametric}
\keyword{regression}

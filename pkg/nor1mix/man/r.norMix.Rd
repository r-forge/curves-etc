\name{r.norMix}
\alias{r.norMix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
r.norMix(obj, x = NULL, xlim = NULL, n = 511, xy.return = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{xy.return}{ ~~Describe \code{xy.return} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of `comp1'}
  \item{comp2 }{Description of `comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj, x = NULL, xlim = NULL, n = 511,
		     xy.return = TRUE, ...)
{
  ## Purpose: Compute r := f / f0; f = normal mixture; f0 = "best" normal approx
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  m <- m.norMix(obj) #-- number of components
  d.o <- dnorMix(obj, x=x, xlim=xlim, n=n)
  dn  <- dnorm(d.o$x, mean = mean.norMix(obj), sd = sqrt(var.norMix(obj)))
  if(xy.return) list(x = d.o$x, y= d.o$y / dn, f0 = dn) else d.o$y / dn
  }
}
\keyword{distribution}
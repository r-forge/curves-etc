#### Examples which use the new feature of more than one 'constraint'.

suppressMessages(library(cobs))

## do *not* show platform info here (as have *.Rout.save), but in 0_pt-ex.R
options(digits = 6)

if(!dev.interactive(orNone=TRUE)) pdf("multi-constr.pdf")

source(system.file("util.R", package = "cobs"))
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> tryCatch.W.E(), showProc.time(), assertError(), relErrV(), ...

Rsq <- function(obj) {
    stopifnot(inherits(obj, "cobs"), is.numeric(res <- obj$resid))
    1 - sum(res^2)/obj$SSy
}
list_ <- function (...) `names<-`(list(...), vapply(sys.call()[-1L], as.character, ""))
is.cobs <- function(x) inherits(x, "cobs")

set.seed(908)
x <- seq(-1,2, len = 50)
f.true <- pnorm(2*x)
y <- f.true + rnorm(50)/10
plot(x,y); lines(x, f.true, col="gray", lwd=2, lty=3)

## constraint on derivative at right end:
(con <- rbind(c(2 , max(x), 0))) # f'(x_n) == 0

## Using 'trace = 3' --> 'trace = 2' inside drqssbc2()

## Regression splines (lambda = 0)
c2   <- cobs(x,y, trace = 3)
c2i  <- cobs(x,y, constraint = c("increase"), trace = 3)
c2c  <- cobs(x,y, constraint = c("concave"), trace = 3)

c2IC <- cobs(x,y, constraint = c("inc", "concave"), trace = 3)
## here, it's the same as just "i":
all.equal(fitted(c2i), fitted(c2IC))

c1   <- cobs(x,y, degree = 1, trace = 3)
c1i  <- cobs(x,y, degree = 1, constraint = c("increase"), trace = 3)
c1c  <- cobs(x,y, degree = 1, constraint = c("concave"), trace = 3)

plot(c1)
lines(predict(c1i), col="forest green")
all.equal(fitted(c1), fitted(c1i), tol = 1e-9)# but not 1e-10

## now gives warning (not error):
c1IC <- cobs(x,y, degree = 1, constraint = c("inc", "concave"), trace = 3)

cp2   <- cobs(x,y,                          pointwise = con, trace = 3)

## Here, warning ".. 'ifl'.. " on *some* platforms (e.g. Windows 32bit) :
r2i <- tryCatch.W.E( cobs(x,y, constraint = "increase", pointwise = con) )
cp2i <- r2i$value
if(doExtras()) print(r2i$warning) # not by default as long as have multi-constr.Rout.save
## when plotting it, we see that it gave a trivial constant!!
cp2c  <- cobs(x,y, constraint = "concave",  pointwise = con, trace = 3)

## now gives warning (not error): but no warning on M1 mac -> IGNORE
## IGNORE_RDIFF_BEGIN
cp2IC <- cobs(x,y, constraint = c("inc", "concave"), pointwise = con, trace = 3)
## IGNORE_RDIFF_END
cp1   <- cobs(x,y, degree = 1,                            pointwise = con, trace = 3)
cp1i  <- cobs(x,y, degree = 1, constraint = "increase",   pointwise = con, trace = 3)
cp1c  <- cobs(x,y, degree = 1, constraint = "concave",    pointwise = con, trace = 3)

cp1IC <- cobs(x,y, degree = 1, constraint = c("inc", "concave"), pointwise = con, trace = 3)

## a named list of all cobs() results above:
cobsL <- mget(Filter(\(nm) is.cobs(.GlobalEnv[[nm]]), ls(patt="c[12p]")),
              envir = .GlobalEnv)

knL <- lapply(cobsL, `[[`, "knots")
str(knL[order(lengths(knL))])

dput(signif(sapply(cobsL, Rsq), digits=8))
Rsqrs <- c(c1  = 0.95079126,  c1c = 0.92974549,
           c1i = 0.95079126, c1IC = 0.92974549,
           c2  = 0.94637437,  c2c = 0.92505977, c2i  = 0.95022829, c2IC = 0.91375404,
           cp1 = 0.9426453,  cp1c = 0.92223149,
           cp1i= 0.9426453, cp1IC = 0.92223149,
           cp2 = 0.94988863, cp2c = 0.91375409, cp2i = 0.93611487,cp2IC = 0.90051964)

all.equal(Rsqrs, sapply(cobsL, Rsq), tolerance=0) #  2.6277e-9 (Lnx F 38)
stopifnot(exprs = {
    all.equal(Rsqrs, sapply(cobsL, Rsq))
    identical(c(5L, 3L, 5L, 3L, 3L, 3L, 4L, 2L,
                5L, 3L, 5L, 3L, 4L, 2L, 4L, 2L),
              unname(lengths(knL)))
})

plot(x,y, main = "cobs(*, degree= 1, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c1),
                 fitted(c1i),
                 fitted(c1c),
                 fitted(cp1),
                 fitted(cp1i),
                 fitted(cp1c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

if(dev.interactive()) x11() # cheap way to get plot in new window, when testing

plot(x,y, main = "cobs(*, degree= 2, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c2),
                 fitted(c2i),
                 fitted(c2c),
                 fitted(cp2),
                 fitted(cp2i),
                 fitted(cp2c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

##--> "increase + pointwise" gives constant which seems plain wrong  <<<< BUG ???

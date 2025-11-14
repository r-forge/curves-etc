library(simest)

(do.pdf <- !dev.interactive(orNone = TRUE))
oop <- options(warnPartialMatchArgs = FALSE, width = 111, digits = 7)

### 1. original from ../man/cvx.lip.reg.Rd  (w/ n = 50) --------------------

mkGauss <-  function(n, sigma, trFn) {
    x <- runif(n, -1,1)
    y <- trFn(x) + sigma * rnorm(n)
    data.frame(x, y)
}

## Approximately geometrically growing n  ("power law"):
n.s <- round(signif(100 * 1.5^(-8:6), 2))
## FIXME:             vvvvvvvvvv  n=1000 etc takes too long !
n.s <- n.s[3 <= n.s & n.s <= 900]
n.s
names(n.s) <- paste0("n=", n.s)

set.seed(47)
d.sqr <- lapply(n.s, mkGauss, sigma = 0.3,  trFn = function(x) x^2)
d.exp <- lapply(n.s, mkGauss, sigma = 0.75, trFn = exp)

## FIXME_2 : the computation is *NOT* interruptable !!!!

prMilli <- function(mat, ...) {
    cat("Times [milli-seconds]:\n")
    print.table(1000 * mat, zero.print = ".", ...)
}


doEst <- function(d, cvxFN, ...) {
    ct <- system.time( fit <- cvxFN(d$x, d$y, ...) # ... e.g. L=L  for cvx.lip.reg()
                      )[1:3]
    print(ct)
    list(sys.tm = ct, fit = fit)
}

doLIP <- function(d, L) doEst(d, cvx.lip.reg, L=L)

lip10.sqr <- lapply(d.sqr, doLIP, L = 10)
lip04.sqr <- lapply(d.sqr, doLIP, L =  4)
lip10.exp <- lapply(d.exp, doLIP, L = 10)
lip04.exp <- lapply(d.exp, doLIP, L =  4)

l10Sqr <- lapply(lip10.sqr, `[[`, "fit")
prMilli(t(sapply(lip10.sqr, `[[`, "sys.tm")))
l04Sqr <- lapply(lip04.sqr, `[[`, "fit")

l10Exp  <-   lapply(lip10.exp, `[[`, "fit")
l04Exp  <-   lapply(lip04.exp, `[[`, "fit")

prMilli(tml10 <- t(sapply(lip10.exp, `[[`, "sys.tm"))) # in milli sec
        tml04 <- t(sapply(lip04.exp, `[[`, "sys.tm"))

aePrt <- function(targ, curr, digits = 3, tol = 0) {
    print(digits=digits,
          with(attributes(all.equal(targ, curr, tolerance = tol, giveErr=TRUE)),
               setNames(err, what)))
}

cat("Timings L=10 vs  L=4:\n"); aePrt(tml10, tml04)

## inspired by sfsmisc::mult.fig() :
## op <- par(mfrow = c(2, 3), las = 1, mgp = c(2, 0.6, 0), mar = 0.1 + c(4, 4, 2, 1))
        ##   oma = c(0, 0, tit.wid, 0), main = NULL, tit.wid = if (is.null(main)) 0 else 1 +
        ## 1.5 * cex.main, cex.main = par("cex.main"), line.main = cex.main -
        ## 1/2, col.main = par("col.main"), font.main = par("font.main"),

if(do.pdf) pdf("cvxLIP_l10_sqr.pdf"); invisible( lapply(l10Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxLIP_l04_sqr.pdf"); invisible( lapply(l10Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxLIP_l10_exp.pdf"); invisible( lapply(l10Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxLIP_l04_exp.pdf"); invisible( lapply(l04Exp, plot) ) ; if(do.pdf) dev.off()

## predict(tmp, newdata = rnorm(10,0,0.1))
set.seed(201)
newd <- unique(sort(0.1 * round(1000*rnorm(250))/1000))
anyDuplicated(newd)
## [1] 18 --- should work, right?

str(p10 <- vapply(l10Sqr, predict, FUN.VALUE=newd, newdata = newd))
p10 # default  7 digits

p04 <- vapply(l04Sqr, predict, FUN.VALUE=newd, newdata = newd)
cat("Predict f=sqrt, L=10 vs L=4:\n"); aePrt(p10, p04)


### 2. original from ../man/cvx.lse.con.reg.Rd  (w/ n = 50)
### 3. and      from ../man/cvx.lse.reg.Rd ---------------------------------------------------

doConreg <- function(d) doEst(d, cvx.lse.con.reg)
doLseReg <- function(d) doEst(d, cvx.lse.reg)

ConR.sqr <- lapply(d.sqr, doConreg) ## warnings:
## 1: In diff(fit)/diff(x) :
##   longer object length is not a multiple of shorter object length

## ==> reason: for n=512 -- two x[]s are _merged_  -> new n = 509 (one w[] = 2) but code
##     ------  uses  diff(fit)/diff(x)  where it should use  diff(fit)/diff(out$x) !!!!! BUG
## BUG  fix must lead to different result structures !! ===> look at smooth.spline() how to do
## =====================================================================
## Also  fastmerge()  should/could profit from smooth.spline() and vice versa ! !
##
ConR.exp <- lapply(d.exp, doConreg)

LseR.sqr <- lapply(d.sqr, doLseReg) # no warnings
LseR.exp <- lapply(d.exp, doLseReg) # no warnings

CRSqr <-  lapply(ConR.sqr, `[[`, "fit")
prMilli(t(sapply(ConR.sqr, `[[`, "sys.tm"))) # in milli sec

CRExp <-  lapply(ConR.exp, `[[`, "fit")
prMilli(t(sapply(ConR.exp, `[[`, "sys.tm")))

if(FALSE) { ## FAILS -- because of the above BUG!
if(do.pdf) pdf("cvxConR_sqr.pdf"); invisible( lapply(CRSqr, plot) ) ; if(do.pdf) dev.off()
}
if(do.pdf) pdf("cvxConR_exp.pdf"); invisible( lapply(CRExp, plot) ) ; if(do.pdf) dev.off()

str(pSq <- vapply(CRSqr, predict, FUN.VALUE=newd, newdata = newd))
pSq # default  7 digits
(pEx <- vapply(CRExp, predict, FUN.VALUE=newd, newdata = newd))

pSq.1 <- vapply(CRSqr, predict, FUN.VALUE=newd, newdata = newd, deriv = 1)
print(pSq.1, digits = 3)

pEx.1 <- vapply(CRExp, predict, FUN.VALUE=newd, newdata = newd, deriv = 1)
print(pEx.1, digits = 3)

## --- LseR

LseRSqr <- lapply(LseR.sqr, `[[`, "fit")
prMilli (t(sapply(LseR.sqr, `[[`, "sys.tm"))) # in milli sec

LseRExp <- lapply(LseR.exp, `[[`, "fit")
prMilli (t(sapply(LseR.exp, `[[`, "sys.tm")))

if(do.pdf) pdf("cvxLseR_sqr.pdf"); invisible( lapply(LseRSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxLseR_exp.pdf"); invisible( lapply(LseRExp, plot) ) ; if(do.pdf) dev.off()

str(pSq <- vapply(LseRSqr, predict, FUN.VALUE=newd, newdata = newd))
pSq # default  7 digits
(pEx <- vapply(LseRExp, predict, FUN.VALUE=newd, newdata = newd))

pSq.1 <- vapply(LseRSqr, predict, FUN.VALUE=newd, newdata = newd, deriv = 1)
print(pSq.1, digits = 3)

pEx.1 <- vapply(LseRExp, predict, FUN.VALUE=newd, newdata = newd, deriv = 1)
print(pEx.1, digits = 3)

### 4. original from ../man/cvx.pen.reg.Rd  (w/ n = 50) ===================================
doPEN <- function(d, lambda, ...) doEst(d, cvx.pen.reg, lambda=lambda, ...)

## here  n=4  produces tons of  NaN  (now gives error  "Number of samples must be at least 5."):
if(!inherits( try(penr <- with(d.sqr[["n=4"]], cvx.pen.reg(x, y, lambda = 0.1))),
             "try-error")) str(penr) ## mostly NaN [C code using `n - 5`]
##=>
d.sqr_ <- d.sqr[n.s >= 5]
d.exp_ <- d.exp[n.s >= 5]

Pen01.sqr <- lapply(d.sqr_, doPEN, lambda = 0.1)
Pen02.sqr <- lapply(d.sqr_, doPEN, lambda = 0.2)
Pen05.sqr <- lapply(d.sqr_, doPEN, lambda = 0.5)
Pen10.sqr <- lapply(d.sqr_, doPEN, lambda = 1.0)

Pen01.exp <- lapply(d.exp_, doPEN, lambda = 0.1)
Pen02.exp <- lapply(d.exp_, doPEN, lambda = 0.2)
Pen05.exp <- lapply(d.exp_, doPEN, lambda = 0.5)

Pen01Sqr <- lapply(Pen01.sqr, `[[`, "fit")
Pen02Sqr <- lapply(Pen02.sqr, `[[`, "fit")
Pen05Sqr <- lapply(Pen05.sqr, `[[`, "fit")
Pen10Sqr <- lapply(Pen10.sqr, `[[`, "fit")
prMilli(t(sapply(Pen01.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(Pen02.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(Pen05.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(Pen10.sqr, `[[`, "sys.tm")))

Pen01Exp <- lapply(Pen01.exp, `[[`, "fit")
Pen02Exp <- lapply(Pen02.exp, `[[`, "fit")
Pen05Exp <- lapply(Pen05.exp, `[[`, "fit")
prMilli(t(sapply(Pen01.exp, `[[`, "sys.tm"))) # in milli sec
prMilli(t(sapply(Pen02.exp, `[[`, "sys.tm"))) # in milli sec
prMilli(t(sapply(Pen05.exp, `[[`, "sys.tm"))) # in milli sec


if(do.pdf) pdf("cvxPen01_sqr.pdf"); invisible( lapply(Pen01Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxPen02_sqr.pdf"); invisible( lapply(Pen02Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxPen05_sqr.pdf"); invisible( lapply(Pen05Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxPen10_sqr.pdf"); invisible( lapply(Pen10Sqr, plot) ) ; if(do.pdf) dev.off()

if(do.pdf) pdf("cvxPen01_exp.pdf"); invisible( lapply(Pen01Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxPen02_exp.pdf"); invisible( lapply(Pen02Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxPen05_exp.pdf"); invisible( lapply(Pen05Exp, plot) ) ; if(do.pdf) dev.off()

dim(prV <- matrix(newd, nrow=length(newd), ncol= ncol(predict(Pen01Sqr[[3]], newdata=1:11))))
str(pSq1 <- vapply(Pen01Sqr, predict, FUN.VALUE=prV, newdata = newd)); pSq1
str(pSq2 <- vapply(Pen02Sqr, predict, FUN.VALUE=prV, newdata = newd)); pSq2
str(pSq5 <- vapply(Pen05Sqr, predict, FUN.VALUE=prV, newdata = newd)); pSq5
str(pSq10<- vapply(Pen10Sqr, predict, FUN.VALUE=prV, newdata = newd)); pSq10
# default  7 digits
(pEx1 <- vapply(Pen01Exp, predict, FUN.VALUE=prV, newdata = newd))
(pEx2 <- vapply(Pen02Exp, predict, FUN.VALUE=prV, newdata = newd))
(pEx5 <- vapply(Pen05Exp, predict, FUN.VALUE=prV, newdata = newd))


### 5. original from ../man/smooth.pen.reg.Rd  (w/ n = 50) ===================================
doSmPen <- function(d, lambda, agcv, ...) doEst(d, smooth.pen.reg, lambda=lambda, agcv=agcv, ...)

SmPen001F.sqr <- lapply(d.sqr, doSmPen, lambda = 0.01, agcv = FALSE)
SmPen001T.sqr <- lapply(d.sqr, doSmPen, lambda = 0.01, agcv = TRUE)
SmPen005F.sqr <- lapply(d.sqr, doSmPen, lambda = 0.05, agcv = FALSE)
SmPen005T.sqr <- lapply(d.sqr, doSmPen, lambda = 0.05, agcv = TRUE)
SmPen02F.sqr <-  lapply(d.sqr, doSmPen, lambda = 0.2,  agcv = FALSE)
SmPen02T.sqr <-  lapply(d.sqr, doSmPen, lambda = 0.2,  agcv = TRUE)

SmPen001F.exp <- lapply(d.exp, doSmPen, lambda = 0.01, agcv = FALSE)
SmPen001T.exp <- lapply(d.exp, doSmPen, lambda = 0.01, agcv = TRUE)
SmPen005F.exp <- lapply(d.exp, doSmPen, lambda = 0.05, agcv = FALSE)
SmPen005T.exp <- lapply(d.exp, doSmPen, lambda = 0.05, agcv = TRUE)
SmPen02F.exp <-  lapply(d.exp, doSmPen, lambda = 0.2,  agcv = FALSE)
SmPen02T.exp <-  lapply(d.exp, doSmPen, lambda = 0.2,  agcv = TRUE)

SmPen001FSqr <- lapply(SmPen001F.sqr, `[[`, "fit")
SmPen001TSqr <- lapply(SmPen001T.sqr, `[[`, "fit")
SmPen001FExp <- lapply(SmPen001F.exp, `[[`, "fit")
SmPen001TExp <- lapply(SmPen001T.exp, `[[`, "fit")

SmPen005FSqr <- lapply(SmPen005F.sqr, `[[`, "fit")
SmPen005TSqr <- lapply(SmPen005T.sqr, `[[`, "fit")
SmPen005FExp <- lapply(SmPen005F.exp, `[[`, "fit")
SmPen005TExp <- lapply(SmPen005T.exp, `[[`, "fit")

SmPen02FSqr <- lapply(SmPen02F.sqr, `[[`, "fit")
SmPen02TSqr <- lapply(SmPen02T.sqr, `[[`, "fit")
SmPen02FExp <- lapply(SmPen02F.exp, `[[`, "fit")
SmPen02TExp <- lapply(SmPen02T.exp, `[[`, "fit")

prMilli(t(sapply(SmPen001F.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen001T.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen001F.exp, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen001T.exp, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen005F.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen005T.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen005F.exp, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen005T.exp, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen02F.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen02T.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen02F.exp, `[[`, "sys.tm")))
prMilli(t(sapply(SmPen02T.exp, `[[`, "sys.tm")))

if(do.pdf) pdf("SmPen001F_sqr.pdf"); invisible( lapply(SmPen001FSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen001T_sqr.pdf"); invisible( lapply(SmPen001TSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen001F_exp.pdf"); invisible( lapply(SmPen001FExp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen001T_exp.pdf"); invisible( lapply(SmPen001TExp, plot) ) ; if(do.pdf) dev.off()

if(do.pdf) pdf("SmPen005F_sqr.pdf"); invisible( lapply(SmPen005FSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen005T_sqr.pdf"); invisible( lapply(SmPen005TSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen005F_exp.pdf"); invisible( lapply(SmPen005FExp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen005T_exp.pdf"); invisible( lapply(SmPen005TExp, plot) ) ; if(do.pdf) dev.off()

if(do.pdf) pdf("SmPen02F_sqr.pdf"); invisible( lapply(SmPen02FSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen02T_sqr.pdf"); invisible( lapply(SmPen02TSqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen02F_exp.pdf"); invisible( lapply(SmPen02FExp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("SmPen02T_exp.pdf"); invisible( lapply(SmPen02TExp, plot) ) ; if(do.pdf) dev.off()

str(pSqF <- vapply(SmPen001FSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqF
str(pSqT <- vapply(SmPen001TSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqT
str(pExF <- vapply(SmPen001FExp, predict, FUN.VALUE=newd, newdata = newd)); pExF
str(pExT <- vapply(SmPen001TExp, predict, FUN.VALUE=newd, newdata = newd)); pExT

str(pSqF <- vapply(SmPen005FSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqF
str(pSqT <- vapply(SmPen005TSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqT
str(pExF <- vapply(SmPen005FExp, predict, FUN.VALUE=newd, newdata = newd)); pExF
str(pExT <- vapply(SmPen005TExp, predict, FUN.VALUE=newd, newdata = newd)); pExT

str(pSqF <- vapply(SmPen02FSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqF
str(pSqT <- vapply(SmPen02TSqr, predict, FUN.VALUE=newd, newdata = newd)); pSqT
str(pExF <- vapply(SmPen02FExp, predict, FUN.VALUE=newd, newdata = newd)); pExF
str(pExT <- vapply(SmPen02TExp, predict, FUN.VALUE=newd, newdata = newd)); pExT


### 6. original from ../man/sim.est.Rd  (w/ n = 50) =========================================
doSIM <- function(d, method, lambda=NULL, L=NULL, nmulti=3, ...) {
    # "BUG":  sim.est(x, y) does *not* work when ncol(d$y) == 1 !! --->
    if(!is.matrix(d$y) || ncol(d$y) < 2)
        stop("[Bug in sim.est():] d$y must be a matrix with at least 2 columns")
    ## for back compatiblity with simest_0.4
    if (method == "cvx.lse.con") { method <- "cvx.lse.reg"; force <- FALSE }
    if (method == "cvx.lse"    ) { method <- "cvx.lse.reg"; force <- TRUE  }
    else force <- FALSE # unused
    doEst(d, sim.est, method=method, force=force, lambda=lambda, L=L, nmulti=nmulti, ...)
}
simMeths <- c("cvx.pen", "cvx.lip", "cvx.lse.con", "cvx.lse", "smooth.pen")

## FIXME some methods need  lambda , other need  L  -- ??

## A bug in sim.est  ==> n = 6 _fails_  with old simest() --- needs ' , drop=FALSE ' in A[.....] :
##  Error in colSums(fit$residuals * fit$deriv * A[, -c(1, 2)]) (from #2) :
##      'x' must be an array of at least two dimensions

try(
SimAll.sqr <- lapply(setNames(,simMeths),
                     function(meth) lapply(d.sqr_, doSIM, method = meth))
)
try(
SimAll.exp <- lapply(setNames(,simMeths),
                     function(meth) lapply(d.exp_, doSIM, method = meth))
)

if(FALSE) { #++++++++++++++++++++++++ NOT YET ++++++++++++++++++++++++++++++++++++++++++++++++++++

SimAllSqr <- lapply(SimAll.sqr, function(sim) lapply(sim, `[[`, "fit"))
SimAllExp <- lapply(SimAll.exp, function(sim) lapply(sim, `[[`, "fit"))

lapply(SimAll.sqr, function(sim) prMilli(t(sapply(sim, `[[`, "sys.tm"))))
lapply(SimAll.exp, function(sim) prMilli(t(sapply(sim, `[[`, "sys.tm"))))


if(do.pdf) pdf("cvxSim01_sqr.pdf"); invisible( lapply(Sim01Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSim02_sqr.pdf"); invisible( lapply(Sim02Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSim05_sqr.pdf"); invisible( lapply(Sim05Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSim10_sqr.pdf"); invisible( lapply(Sim10Sqr, plot) ) ; if(do.pdf) dev.off()

if(do.pdf) pdf("cvxSim01_exp.pdf"); invisible( lapply(Sim01Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSim02_exp.pdf"); invisible( lapply(Sim02Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSim05_exp.pdf"); invisible( lapply(Sim05Exp, plot) ) ; if(do.pdf) dev.off()

## predict
str(pSq1 <- vapply(Sim01Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq1
str(pSq2 <- vapply(Sim02Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq2
str(pSq5 <- vapply(Sim05Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq5
str(pSq10<- vapply(Sim10Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq10

(pEx1 <- vapply(Sim01Exp, predict, FUN.VALUE=newd, newdata = newd))
(pEx2 <- vapply(Sim02Exp, predict, FUN.VALUE=newd, newdata = newd))
(pEx5 <- vapply(Sim03Exp, predict, FUN.VALUE=newd, newdata = newd))

}### NOT YET ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### 7. original from ../man/simestgcv.Rd  (w/ n = 50) ==========================================
doSimCV <- function(d, lambda, ...) doEst(d, simestgcv, lambda=lambda, ...)

try(
SimCV01.sqr <- lapply(d.sqr, doSimCV, lambda = 0.1) # error invalid 'xmax'
)

if(FALSE) { #++++++++++++++++++++++++ NOT YET ++++++++++++++++++++++++++++++++++++++++++++++++++++

SimCV02.sqr <- lapply(d.sqr, doSimCV, lambda = 0.2)
SimCV05.sqr <- lapply(d.sqr, doSimCV, lambda = 0.5)
SimCV10.sqr <- lapply(d.sqr, doSimCV, lambda = 1.0)

SimCV01.exp <- lapply(d.exp, doSimCV, lambda = 0.1)
SimCV02.exp <- lapply(d.exp, doSimCV, lambda = 0.2)
SimCV05.exp <- lapply(d.exp, doSimCV, lambda = 0.5)

SimCV01Sqr <- lapply(SimCV01.sqr, `[[`, "fit")
SimCV02Sqr <- lapply(SimCV02.sqr, `[[`, "fit")
SimCV05Sqr <- lapply(SimCV05.sqr, `[[`, "fit")
SimCV10Sqr <- lapply(SimCV10.sqr, `[[`, "fit")
prMilli(t(sapply(SimCV01.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SimCV02.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SimCV05.sqr, `[[`, "sys.tm")))
prMilli(t(sapply(SimCV10.sqr, `[[`, "sys.tm")))

SimCV01Exp <- lapply(SimCV01.exp, `[[`, "fit")
SimCV02Exp <- lapply(SimCV02.exp, `[[`, "fit")
SimCV05Exp <- lapply(SimCV05.exp, `[[`, "fit")
prMilli(t(sapply(SimCV01.exp, `[[`, "sys.tm"))) # in milli sec
prMilli(t(sapply(SimCV02.exp, `[[`, "sys.tm"))) # in milli sec
prMilli(t(sapply(SimCV05.exp, `[[`, "sys.tm"))) # in milli sec


## if(FALSE) { ## FAILS -- because of the above BUG!
if(do.pdf) pdf("cvxSimCV01_sqr.pdf"); invisible( lapply(SimCV01Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSimCV02_sqr.pdf"); invisible( lapply(SimCV02Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSimCV05_sqr.pdf"); invisible( lapply(SimCV05Sqr, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSimCV10_sqr.pdf"); invisible( lapply(SimCV10Sqr, plot) ) ; if(do.pdf) dev.off()
## }
if(do.pdf) pdf("cvxSimCV01_exp.pdf"); invisible( lapply(SimCV01Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSimCV02_exp.pdf"); invisible( lapply(SimCV02Exp, plot) ) ; if(do.pdf) dev.off()
if(do.pdf) pdf("cvxSimCV05_exp.pdf"); invisible( lapply(SimCV05Exp, plot) ) ; if(do.pdf) dev.off()

str(pSq1 <- vapply(SimCV01Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq1
str(pSq2 <- vapply(SimCV02Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq2
str(pSq5 <- vapply(SimCV05Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq5
str(pSq10<- vapply(SimCV10Sqr, predict, FUN.VALUE=newd, newdata = newd)); pSq10
# default  7 digits
(pEx1 <- vapply(SimCV01Exp, predict, FUN.VALUE=newd, newdata = newd))
(pEx2 <- vapply(SimCV02Exp, predict, FUN.VALUE=newd, newdata = newd))
(pEx5 <- vapply(SimCV03Exp, predict, FUN.VALUE=newd, newdata = newd))

}### NOT YET ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






options(oop) ## revert options() at end

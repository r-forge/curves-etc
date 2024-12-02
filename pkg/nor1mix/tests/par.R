library(nor1mix)

is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()[["sysname"]] == "Darwin" # M1mac etc  are less accurate

## Check  nM2par(), par2norMix() and llnorMix() :
op <- options(keep.source = FALSE) # to keep this small
chk <- quote({
        all.equal(pp,   nM2par(nm),                   tol= tol1)
        all.equal(pp.l, nM2par(nm.l, trafo= "logit"), tol= tol1)
        all.equal(pp.c, nM2par(nm.c, trafo=  "clr1"), tol= tol1)
        all.equal(obj, nm,   check.attributes=FALSE, tol= tol2)
        all.equal(obj, nm.l, check.attributes=FALSE, tol= tol2)
        all.equal(obj, nm.c, check.attributes=FALSE, tol= tol2)
        ## xx
        all.equal(llnorMix(pp  , xx),                 logLik.x, tol = tol1)
        all.equal(llnorMix(pp.l, xx, trafo= "logit"), logLik.x, tol = tol1)
        all.equal(llnorMix(pp.c, xx, trafo=  "clr1"), logLik.x, tol = tol1)
}); options(op)
chkTrue <- do.call(expression, as.list(chk)[-1L])

tol1 <- (if(isMac) 4 else 1) * 1e-15
tol2 <- (if(isMac) 8 else 4) * 1e-15
nms <- paste("MW.nm", 1:16, sep="")
for(n in nms) {
    cat(n,":")
    obj <- get(n, envir = as.environment("package:nor1mix"))
    xx <- rnorMix(1000, obj)
    logLik.x <- sum(dnorMix(xx, obj, log = TRUE))
    pp   <- nM2par(obj) # use "current default"
    pp.l <- nM2par(obj, trafo = "logit")
    pp.c <- nM2par(obj, trafo = "clr1")
    nm   <- par2norMix(pp) # (current) default
    nm.l <- par2norMix(pp.l, trafo= "logit")
    nm.c <- par2norMix(pp.c, trafo=  "clr1")
    ## first show (tol* = 0), then check with the tolerances {tol1, tol2}
    local({ tol1 <- tol2 <- 0
        cat("With  tol1 = tol2 = 0 :\n")
        for(ee in chkTrue)
            cat(deparse1(ee), ": ", eval(ee), "\n")
    })
    stopifnot(exprObject = chkTrue)
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
}

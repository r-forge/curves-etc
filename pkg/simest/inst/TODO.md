1. A bug: whenever `fastmerge()` (or `smooth.spline()`) is used, the data of
  too close x-values may be merged.
  Currently, cvx.lse.conreg() in `diff(fit)/diff(x)` does not take this
  into account ... but at least  `conreg()` itself, does so
  (using an older version of "merge()", that smooth.spline() used
   originally: `xx <- signif(x, 6))`

2. [related to `1.`]:  `fastmerge()` does "something similar" as the
   initial part of  `smooth.spline()`:  `~/R/D/r-devel/R/src/library/stats/R/smspline.R`
   ---> MM thinks they __should learn from each other__

3. Numerical stability:   `diff(y)/diff(x)` is really "simplistic"
  Considerably better is what  `stats:::splinefun(...., method = "monoH.FC")` uses,
  viz. _"central differences"_ in the middle: the  `(Sx[-1] + Sx[-n1])/2` :

```{r}
  nx <- length(x)
  n1 <- nx - 1L  # stopifnot(nx > 0, !anyNA(x))

  dy <- y[-1L] - y[-nx] # = diff(y), slightly faster
  dx <- x[-1L] - x[-nx] # = diff(x), slightly faster
  Sx <- dy/dx
  m <- c(Sx[1L], (Sx[-1L] + Sx[-n1])/2, Sx[n1])
```

4. `sim.est(x, y)`  needs  `ncol(x) >= 2`  (such that a norm-1 "beta" makes sense!)
   ==> Improve such that for ncol(x) == 1,  beta is *not* estimated (but
   beta :== 1 constant).

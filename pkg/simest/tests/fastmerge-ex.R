library(simest)
## relevant example --- used in ../man/fastmerge.Rd
n <- 47 ; j <- 0L
while(j < 1e4) {
    set.seed(j)
    x <- signif(runif(n, -1,1), 5)
    y <- signif(runif(n, -1,1), 5)
    DataMat <- cbind(x, y)
    fmLst <- fastmerge(DataMat)
    if((nw <- length(w <- fmLst$w)) < n) {
        cat("found (x,y) -- j =", j," nw =",nw,"\n")
        if((nX <- n - nw) > 1) {
            cat("** n - |w| =", nX, "--- table(w):\n")
            print(table(w))
        }
    }
    j <- j + 1L
}
## gives a lot ... and then
## found (x,y) -- j = 2657  nw = 44 
## ** n - |w| = 3 --- table(w):
## w
##  1  2 
## 41  3 

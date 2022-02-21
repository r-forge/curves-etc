library("plugdensity")

options(digits = 6)

data(faithful)
(pd.geys <- plugin.density(faithful$waiting))
pd.geys$y * 1e4

data(co2)
(pd.co2 <- plugin.density(co2))

## Add one outlier to co2, move it to Inf --- bw.EH() is quite *robust*
xOs <- c(max(co2), 370+ 2^c(0:8, 16*(1:6)), Inf)
bwOs <- vapply(xOs, function(xo) bw.EH( c(co2, xo)), 0.1)
bwOm <- vapply(xOs, function(xo) bw.EH(-c(co2, xo)), 0.1)
cbind(xOs, bwOs, D = bwOs - bwOs[length(xOs)], bwOm)
plot(bwOs ~ xOs, log="x", type="b")
plot (bwOs ~ xOs, log="x", type="b", subset = xOs < 4e4)
lines(bwOm ~ xOs, col=2, type="b")

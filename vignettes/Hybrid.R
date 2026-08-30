## -----------------------------------------------------------------------------
library(ramp.falciparum)
library(ramp.func)
library(deSolve)
library(knitr)

## -----------------------------------------------------------------------------
# devtools::load_all()

## ----fig.height=4, fig.width=7------------------------------------------------
m = seq(0.1, 24)
N = pmax(24, 4*m)
plot(m, ppois(N,m), type = "l", ylim = c(1-1e-07,1+1e-07))

## ----fig.height=4, fig.width=7, echo=F----------------------------------------
FoI_a = make_F_a(5/365, age_par = makepar_F_type2(), season_par=makepar_F_sin())

a2years  = 0:730

plot(a2years, FoI_a(a2years), type = "l", 
     xlab = "a - age (in days)", ylab = expression(FoI[d](a)))

## -----------------------------------------------------------------------------
solveMMinfty(3/365, FoI_a) -> out

## ----fig.height=4, fig.width=7, echo=F----------------------------------------
MoIDistPlot(out$moi, 500)

## ----fig.height=4, fig.width=7------------------------------------------------
solve_dm(3/365, FoI_a) -> mt 
plot(mt$time, mt$m, type = "l", 
     xlab = "Age (in Days)", ylab = expression(m[d](a)))

## ----compute 3 ways-----------------------------------------------------------
aa = 0:1825
v1 = meanMoI(aa, FoI_a, hhat=5/365)
v2 = solveMMinfty(5/365, FoI_a, Amax=1825)$m
v3 = solve_dm(5/365, FoI_a, Amax=1825)[,2]

## ----fig.height=7, fig.width=8------------------------------------------------
par(mfrow = c(2,2))

plot(aa, v1, type = "l", col = "red", lwd=4, 
     xlab = "Age (in Days)", 
     ylab = expression(m[d](a)), 
     main = "meanMoI")

plot(aa, v2, type = "l", lwd=2, 
     xlab = "Age (in Days)", 
     ylab = expression(m[d](a)), 
     main = "solveMMinfty")
lines(aa, v2, type = "l", col = "yellow", lwd=2, lty=2)

plot(aa, v3, type = "l", col = "darkblue", lwd=1, 
     xlab = "Age (in Days)", 
     ylab = expression(m[d](a)), 
     main = "solve_dm")

plot(aa, v1, type = "l", col = "red", lwd=4, 
     xlab = "Age (in Days)", 
     ylab = expression(m[d](a)), 
     main = "All Three")
lines(aa, v2, type = "l", col = "yellow", lwd = 2)
lines(aa, v3, type = "l", col = "darkblue", lwd = 2, lty =2)

## ----fig.height=4, fig.width=6, echo=F----------------------------------------
tm = 0:1825
Xt= solve_dpda(5/365, FoI_a, Amax=5*365)
tPR = truePR(tm, FoI_a, hhat=5/365)
plot(tm, tPR, type = "l", ylab = expression(p(a)), xlab = "Age (Days)", lwd=5, col = "red")
lines(Xt$time, 1-exp(-Xt$m), col = "white", lwd=2)
lines(Xt$time, Xt$p, col = "blue", lwd=2, lty=2)

## ----fig.height=5, fig.width=7------------------------------------------------
par(mar = c(7, 4, 2, 2))
solve_dAoI(5/365, FoI_a) -> mt 
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(list(x, sqrt(x[paste("[2]")]), sqrt(x[paste("[3]")], 3))), ylim = range(mt$xn3^(1/3))) 
lines(mt$time, mt$xn2^(1/2), col = "purple")
lines(mt$time, mt$xn3^(1/3), col = "darkgreen") 

## -----------------------------------------------------------------------------
aa = seq(5, 5*365, by = 5) 
moment1 = momentAoI(aa, FoI_a)
moment2 = momentAoI(aa, FoI_a, n=2)
moment3 = momentAoI(aa, FoI_a, n=3)

## -----------------------------------------------------------------------------
solve_dAoI(5/365, FoI_a, Amax = 5*365, dt=5) -> mt 

## ----fig.height=7, fig.width=7------------------------------------------------
par(mfrow = c(2,2), mar = c(7, 4, 2, 2))
# Top Left
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(x), lwd=3) 
lines(aa, moment1, col = "yellow", lwd=2, lty=2) 
# Top Right 
plot(mt$time, mt$xn2, type = "l", xlab = "Age (in Days)", ylab = expression(x[paste("[2]")]), lwd=3, col = "purple") 
lines(aa, moment2, col = "yellow", lwd=2, lty=2)
plot(mt$time, mt$xn3, type = "l", lwd=3, col = "darkgreen",
     xlab = "Age (in Days)", ylab = expression(x[paste("[3]")]))
lines(aa, moment3, col = "yellow", lwd=2, lty=2)
par(mar = c(7, 4, 2, 2))
solve_dAoI(5/365, FoI_a) -> mt 
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(list(x, sqrt(x[paste("[2]")]), sqrt(x[paste("[3]")], 3))), lwd=3, ylim = range(mt$xn3^(1/3))) 
lines(mt$time, mt$xn2^(1/2), col = "purple", lwd=3)
lines(mt$time, mt$xn3^(1/3), col = "darkgreen", lwd=3) 
lines(aa, moment1, col = "yellow", lwd=2, lty=2)
lines(aa, moment2^(1/2), col = "yellow", lwd=2, lty=2)
lines(aa, moment3^(1/3), col = "yellow", lwd=2, lty=2)

## -----------------------------------------------------------------------------
solve_dAoYda(5/365, FoI_a, Amax=5*365, dt=5, n=9) -> mt

## -----------------------------------------------------------------------------
moment1y = momentAoY(aa, FoI_a, hhat=5/365)

## ----fig.height=4, fig.width=6, echo=F----------------------------------------
plot(aa, moment1y, type = "l", lwd=2, xlab = "Age (in Days)")
with(mt, lines(time, y, col = "darkred", lwd=2))

## ----fig.height=4, fig.width=6, echo=F, eval=F--------------------------------
# plot(moment1y-mt$y[-1], type = "l", xlab = "Age (in Days)", ylab = "Errors")


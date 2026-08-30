## -----------------------------------------------------------------------------
library(ramp.falciparum)
library(ramp.func)
library(deSolve)
library(knitr)

## -----------------------------------------------------------------------------
# devtools::load_all()

## ----echo=F-------------------------------------------------------------------
aa = 0:1095
FoI_a = make_F_a(5/365, age_par = makepar_F_type2(), season_par=makepar_F_sin())

## ----fig.height=4, fig.width=7, echo=F, eval=F--------------------------------
# plot(aa, FoI_a(aa), type = "l",
#      xlab = "a - age (in days)", ylab = expression(FoI[d](a)))

## -----------------------------------------------------------------------------
MMinf <- solveMMinfty(5/365, FoI_a, Amax=1095)


## ----fig.height=4, fig.width=7------------------------------------------------
with(MMinf, plot(time, m, type = "l", ylab = expression(m[d](a)), xlab = "a - cohort age (in days)"))

## -----------------------------------------------------------------------------
hybrid = solve_dm(5/365, FoI_a, Amax=1095)

## ----fig.height=4, fig.width=7------------------------------------------------
with(hybrid, plot(time, m, type = "l", ylab = expression(m[d](a)), xlab = "a - cohort age (in days)"))

## -----------------------------------------------------------------------------
moi = meanMoI(aa, FoI_a, hhat=5/365)

## ----fig.height=4, fig.width=7------------------------------------------------
plot(aa, moi, type = "l", ylab = expression(m[d](a)), xlab = "a - cohort age (in days)")

## -----------------------------------------------------------------------------
c(mean(abs(moi - hybrid$m)) < 1e-9,
mean(abs(MMinf$m- hybrid$m)) < 1e-10,
max(abs(moi - hybrid$m)) < 1e-8,
max(abs(MMinf$m- hybrid$m))< 1e-10)


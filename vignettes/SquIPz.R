## -----------------------------------------------------------------------------
library(ramp.falciparum)
library(ramp.func)

## -----------------------------------------------------------------------------
# devtools::load_all("~/git/ramp.func")
# devtools::load_all()

## ----fig.width=7, fig.height=4------------------------------------------------
F_z = make_NBtrunc(1, 2)
plot(F_z(1:10), type = "h")

## ----fig.width=7--------------------------------------------------------------
F_a = make_F_a(1)
q_out <- solve_SquIPz(8/365, F_a, F_z, sigma=2/365, xi=0, Amax=3*365, N=400)

## ----fig.width=7, fig.height=4------------------------------------------------
with(q_out, plot(age, m_1, type = "l", ylab = expression(m[1]), main = "First Moment of the MoI"))
with(q_out, lines(age, m1, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4------------------------------------------------
ylm = 1e-13
with(q_out, plot(age, m_1-m1, type = "l", ylab = "Errors", main = "First Moment: Numerical Errors", ylim = c(-ylm, ylm)))

## ----fig.width=7, fig.height=4------------------------------------------------
with(q_out, plot(age/365, m_2, type = "l", main = "Second Moment of the MoI",
                 ylab = expression(m[2]), 
                 xlab = "Age", ylim = range(m_2, m2)))
with(q_out, lines(age/365, m2, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4------------------------------------------------
ylm = 1e-10/5
with(q_out, plot(age, m_2-m2, type = "l", 
                 ylab = "Errors", 
                 main = "Second Moment: Numerical Errors", 
                 ylim = c(-ylm, ylm)))

## ----fig.width=7, fig.height=4------------------------------------------------
q1_out <- solve_SquIPz(8/365, F_a, F_z, sigma=2/365, xi=0, Amax=3*365, N=4)
with(q1_out, plot(age, m_1, type = "l", 
                  ylab = expression(m[1]), 
                  main = "First Moment of the MoI, N=4"))
with(q1_out, lines(age, m1, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4------------------------------------------------
clrs = viridisLite::turbo(7)
set.seed(234)
Sa = makepar_F_type2()
Sp = makepar_F_sin()
Tp = makepar_F_spline(seq(0, 3650, length.out=11), 1+runif(11, -1, 1), X=2)
Kp = makepar_F_sharkbite(D=730, L=365)
F_t <- make_ts_function(scale = 0.05, season_par=Sp, trend_par=Tp, shock_par=Kp)
F_a <- make_F_a(avg = 3/365, age_par=Sa, season_par=Sp, trend_par=Tp, shock_par=Kp)

tt <- seq(0, 3650, by=5)
aa <-seq(0, 365*5, by =5) 
plot(tt, 0.05*F_t(tt), type = "l")

## ----fig.width=7, fig.height=4------------------------------------------------
plot(tt, 0.1*F_t(tt), type = "l", lwd=2, ylim = c(0,0.04))
lines(aa, F_a(aa), type = "l", col = clrs[2])
lines(aa+365, F_a(aa, 365), type = "l", col = clrs[3])
lines(aa+730, F_a(aa, 730), type = "l", col = clrs[4])
lines(aa+1095, F_a(aa, 1095), type = "l", col = clrs[5])
lines(aa+1460, F_a(aa, 1460), type = "l", col = clrs[6])

## ----fig.width=7, fig.height=4------------------------------------------------
plot(aa, F_a(aa), type = "n", lwd=2, ylim = c(0,0.04))
lines(aa, F_a(aa), type = "l", col = clrs[2])
lines(aa, F_a(aa, 365), type = "l", col = clrs[3])
lines(aa, F_a(aa, 730), type = "l", col = clrs[4])
lines(aa, F_a(aa, 1095), type = "l", col = clrs[5])
lines(aa, F_a(aa, 1460), type = "l", col = clrs[6])


## ----fig.width=7, fig.height=4------------------------------------------------
plot(aa, cumsum(F_a(aa)), ylab = "Age", xlab = "Cumulative Exposure", type = "n", ylim = c(0,3))
lines(aa, cumsum(F_a(aa)), type = "l", col = clrs[2])
lines(aa, cumsum(F_a(aa, 365)), type = "l", col = clrs[3])
lines(aa, cumsum(F_a(aa, 730)), type = "l", col = clrs[4])
lines(aa, cumsum(F_a(aa, 1095)), type = "l", col = clrs[5])
lines(aa, cumsum(F_a(aa, 1460)), type = "l", col = clrs[6])

## ----fig.width=7, fig.height=4------------------------------------------------
q3_out1 <- solve_SquIPz(5/365, F_a, F_z, sigma=2/365, xi=0, Amax=5*365)
q3_out2 <- solve_SquIPz(5/365, bday=365, F_a, F_z, sigma=2/365, xi=0, Amax=5*365)
q3_out3 <- solve_SquIPz(5/365, bday=730, F_a, F_z, sigma=2/365, xi=0, Amax=5*365)
q3_out4 <- solve_SquIPz(5/365, bday=365*3, F_a, F_z, sigma=2/365, xi=0, Amax=5*365)
q3_out5 <- solve_SquIPz(5/365, bday=365*4, F_a, F_z, sigma=2/365, xi=0, Amax=5*365)
with(q3_out1, plot(age, m_1/x, col = clrs[2], ylab = "mean MoI", type ="l", ylim = c(1, 2)))
with(q3_out2, lines(age, m_1/x, col = clrs[3]))
with(q3_out3, lines(age, m_1/x, col = clrs[4]))
with(q3_out4, lines(age, m_1/x, col = clrs[5]))
with(q3_out5, lines(age, m_1/x, col = clrs[6]))


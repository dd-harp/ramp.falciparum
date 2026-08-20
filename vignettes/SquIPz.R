## -----------------------------------------------------------------------------
library(ramp.falciparum)
devtools::load_all()

## -----------------------------------------------------------------------------
F_z = function(x, mu=1, size=2){return(dnbinom(x,mu=mu, size=size)/(1-dnbinom(0,mu=mu,size=size)))}
plot(F_z(1:10), type = "h") 

## ----fig.width=7--------------------------------------------------------------
foiP1 = list(hbar = 1, agePar = par_flatAge(), 
             seasonPar = par_flatSeason(), 
             trendPar = par_flatTrend())

q_out <- solve_SquIPz(8/365, foiP1, F_z, sigma=2/365, xi=0, Amax=3*365, N=400)

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
q1_out <- solve_SquIPz(8/365, foiP1, F_z, sigma=2/365, xi=0, Amax=3*365, N=4)
with(q1_out, plot(age, m_1, type = "l", 
                  ylab = expression(m[1]), 
                  main = "First Moment of the MoI, N=4"))
with(q1_out, lines(age, m1, lty=2, col = "yellow"))


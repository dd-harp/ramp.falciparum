## -----------------------------------------------------------------------------
library(ramp.falciparum)

## ----fig.width=7, eval=F------------------------------------------------------
# foiP1 = list(hbar = 1, agePar = par_flatAge(), seasonPar = par_flatSeason(), trendPar = par_flatTrend())
# q_out <- solve_SquIP(8/365, sigma=2/365, xi=0, foiP1, Amax=5*365)

## ----fig.width=7, fig.height=4, eval=F----------------------------------------
# with(q_out, plot(age, m, type = "l"))
# with(q_out, lines(age, mm, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4,eval=F-----------------------------------------
# with(q_out, plot(age, m2, type = "l", ylim = range(m2, mm2)))
# with(q_out, lines(age, mm2, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4, eval=F----------------------------------------
# m_out <- solveSIPm(8/365, sigma=2/365, xi=0, foiP1, Amax=5*365)
# with(m_out, plot(age, S, ylim = c(0,1000), type = "l"))
# with(m_out, lines(age, I))
# with(m_out, lines(age, P))
# with(q_out, lines(age, S, col = "yellow", lty=2))
# with(q_out, lines(age, I, col = "yellow", lty=2))
# with(q_out, lines(age, P, col = "yellow", lty=2))

## ----fig.width=7, fig.height=4, eval=F----------------------------------------
# with(q_out, plot(age, Fr, type = "l", ylim = c(0, 0.005)))
# with(m_out, lines(age, c(.005, Fr), col = "red"))

## ----fig.width=7, fig.height=4, eval=F----------------------------------------
# with(m_out, plot(age, m, type = "l"))
# with(q_out, lines(age, m, lty=2, col = "yellow"))

## ----fig.width=7, fig.height=4, eval=F----------------------------------------
# with(m_out, plot(age, m2, type = "l"))
# with(q_out, lines(age, m2, lty=2, col = "yellow"))


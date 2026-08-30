# SIPm

This describes a queuing model `SquIP`

> see `SIP.R`

``` r

library(ramp.falciparum)
library(ramp.func)
```

## SquIP: The Master Equations

`SquIP` is an infinite system of differential equations describing the
dynamics of the distribution of the MoI in a cohort of humans as it
ages.

The system has three classes of variable:

- $`S`$ describes the density of susceptible individuals

- $`P`$ describes the density of chemo-protected individuals

- $`I_i`$ describes the density of individuals infected with exactly
  $`i`$ parasite types

The model has several parameters that could vary with respect to age:

- $`\mu`$ is the death rate of the cohort with respect to age

- $`h`$ is the force of infection (FoI), the rate that infections occur

- $`r`$ is the rate that each clone clears, and we assume that each
  parasite clone clears independently of the others at a constant rate,
  regardless of the age of the clone

- $`\rho`$ is the fraction of incident infections that are treated and
  cured

- $`\xi`$ is the rate that people take drugs and cure infections without
  regard to whether they are infected

- $`\sigma`$ is the rate that people take drugs and cure infections
  because they are infected, so that anyone who is infected would cure
  infections at the rate $`\sigma + \xi`$

- $`\eta`$ is the rate that chemoprotection is lost after curing an
  infection

### Dynamics

The dynamics of susceptible and chemoprotected individuals are:

``` math
\begin{equation}
    \begin{array}{rl}
    \frac{dP}{da} &= - \eta P + (h \rho + \xi) \left( S + I \right) + \sigma I  -\mu P\\ 
    \frac{dS}{da} &= \eta P + r I_1  - (h + \mu + \xi) S
    \end{array}
\end{equation}
```

We track the number infected with $`i`$ types, $`I_i.`$

``` math
\begin{equation}
    \begin{array}{rl}
    \frac{dI_1}{da} &=  2r I_2  + h (1-\rho) S - (r + h + \sigma + \xi + \mu ) I_1 \\ 
    \frac{dI_2}{da} &=  3r I_3 + h (1-\rho) I_1 - (2r + h + \sigma + \xi + \mu) I_2 \\ 
    \frac{dI_i}{da} &=  ir I_i + h (1-\rho) I_{i-1} - (ir + h + \sigma + \xi + \mu) I_i \\ 
    \end{array}
\end{equation}
```

The implementation truncates the infinite system, choosing a largest
MoI, $`n`$:

``` math
\begin{equation}
\frac{dI_n}{da} = h (1-\rho) I_{n-1} - (n r + h \rho + \sigma + \xi+ \mu) I_n
\end{equation}
```

In this implementation, the function `solve_SquIP` solves the equations
and then it computes the first two moments of the distribution:

``` math
m_1 = \frac 1 H \sum_i i I_i
```

and

``` math
m_2 = \frac 1 H \sum_i i^2 I_i
```

In `SIPm` (below), we derive equations to compute $`dm/da`$ and
$`d m_2 /da.`$

denote the mean MoI in the population. We can show that

``` math
 
\frac{\textstyle{dm}}{\textstyle{da}} = h(1-\rho)\left(1-\frac{P}{H} \right) - \left(r + \rho h + \sigma + \xi  \right)  m
```
We also compute the 2nd moment:

``` math
 
\frac{\textstyle{dm_2}}{\textstyle{da}} = h(1-\rho)\left(1-\frac{P}{H} + 2m\right) - \left(2r + \rho h + \sigma + \xi  \right)  m_2 + r m
```

The function `solveSIPqueue`:

- calls `d_SIPqueue_da` which computes all the derivaves for $`S`$,
  $`I`$, $`P`$, and $`m`$

- computes $`m`$ from $`I`$

``` r

foiP1 = list(hbar = 1, agePar = par_flatAge(), seasonPar = par_flatSeason(), trendPar = par_flatTrend())
q_out <- solve_SquIP(8/365, sigma=2/365, xi=0, foiP1, Amax=5*365)
```

For verification, the variable $`m`$ is computed two ways We note that
the variable $`m`$ matches the MoI computed from $`I`$:

``` r

with(q_out, plot(age, m, type = "l"))
with(q_out, lines(age, mm, lty=2, col = "yellow"))
```

``` r

with(q_out, plot(age, m2, type = "l", ylim = range(m2, mm2)))
with(q_out, lines(age, mm2, lty=2, col = "yellow"))
```

## Approximating Equations

``` math
\begin{equation}
    \begin{array}{rl}
    \frac{dP}{da} &= - \eta P + (h \rho + \xi) \left( S + I \right) + \sigma I  -\mu P\\ 
    \frac{dS}{da} &= \eta P + F_r(m, m_2) I  - (h + \mu + \xi) S \\
    \frac{dI}{da} &= h(1-\rho) S - F_r((m, m_2)) I  - (h \rho + \sigma + \xi) S - \mu S
    \end{array}
\end{equation}
```

We compute $`m_1`$ and $`m_2`$ as before, and we let

``` math
F_r = r \frac{\mbox{dnbinom}(1, m, \mbox{size} = \phi)}{1-\mbox{dnbinom}(0, m, \mbox{size} =\phi))}\left(1-\frac P H \right)
```

where
``` math
\phi = \frac{m^2}{m_2-m^2-m}
```

``` r

m_out <- solveSIPm(8/365, sigma=2/365, xi=0, foiP1, Amax=5*365)
with(m_out, plot(age, S, ylim = c(0,1000), type = "l"))
with(m_out, lines(age, I))
with(m_out, lines(age, P))
with(q_out, lines(age, S, col = "yellow", lty=2))
with(q_out, lines(age, I, col = "yellow", lty=2))
with(q_out, lines(age, P, col = "yellow", lty=2))
```

``` r

with(q_out, plot(age, Fr, type = "l", ylim = c(0, 0.005)))
with(m_out, lines(age, c(.005, Fr), col = "red"))
```

``` r

with(m_out, plot(age, m, type = "l"))
with(q_out, lines(age, m, lty=2, col = "yellow"))
```

``` r

with(m_out, plot(age, m2, type = "l"))
with(q_out, lines(age, m2, lty=2, col = "yellow"))
```

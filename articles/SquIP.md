# SquIP

------------------------------------------------------------------------

**SquIP** is a dynamical system describing the dynamics of the
multiplicity of infection (MoI) in a cohort of humans as it ages.
In-Line documentation for **SquIP** is accessed as:

    help(SquIP)

**SquIP** was developed by extending the master equations for the
queuing process model $`M/M/\infty,`$ by adding terms that describe
treating and curing infections and a chemoprotected class.

**Implementation** — The file `SquIP.R` defines three functions:

    dSquIP 
    solve_SquIP
    parse_SquIP

This numerical implementation truncates the infinite system of
differential equations. To verify the code and check the accuracy of the
truncated system, we also derive and compute hybrid variables describing
the first few moments of the infinite system.

**Related**

- [SquIPz](https://dd-harp.github.io/ramp.falciparum/articles/SquIPz.md)
  extends **SquIP** with a Tweedie process for heterogeneous exposure.

- SIPm

------------------------------------------------------------------------

## The Model

### Variables

The independent variable is age, $`a.`$

The dependent variables are:

- $`S`$ is the expected number uninfected and susceptible to infection

- $`P`$ is the expected number uninfected and chemo-protected

- $`I_i`$ is the expected number infected with $`i`$ clones

We let $`H`$ denote the expected population density of the cohort as it
ages:
``` math
H = S+P+\sum_i I_i.
```

Similarly, we introduce the variable $`I`$:

``` math
I = \sum_i I_i
```

### Initial Conditions

In the model, we assume that everyone is born susceptible. $`H`$ is
passed as a parameter, and the initial conditions are set to

- $`S(0) = H`$

- $`P(0) = 0`$

- $`I_i(0) = 0`$

### Mortality

The model assumes that all individuals die as they age at the rate
$`\mu`$, such that

``` math
\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H
```

### Exposure and Infection

We let $`h`$ denote the force of infection (FoI). This model, **SquIP,**
assumes that each incident infection increases the MoI by exactly one.
The model **SquIPz** implements exposure as a Tweedie process.

The model assumes that each clone clears independently of the others at
rate $`r,`$ so for those infected with $`i`$ clones, the MoI goes down
by one at the rate $`ri.`$

### Treatment and Chemprotection

These equations assume that individuals are treated and cured for
several different reasons:

- Incident infections cause disease and a fraction $`\rho`$ gets treated

- Everyone in the population takes drugs at a background rate $`\xi`$

- In addition to other modes of treatment, infected individuals get
  treated at the higher rate $`\sigma`$

After being treated and cured, individuals enter the chemoprotected
class $`P.`$ Chemoprotection is lost at the rate $`\eta,`$ and they
enter the susceptible class $`S.`$

### Differential Equations:

Changes in the susceptible population are:

``` math
\frac{\textstyle{dS}}{\textstyle{da}}  = \eta P + r I_1 - (h + \xi + \mu) S
```

Changes in the chemoprotected population are:

``` math
\frac{\textstyle{dP}}{\textstyle{da}}  =  \left(h \rho + \xi \right) \left(H-P \right) + \sigma \sum_i I_i   - \left(\eta + \mu \right) P
```

Changes in the infected population are given by an infinite system of
equations. The first one is a special case:

``` math
\frac{\textstyle{dI_1}}{\textstyle{da}}  =  h(1-\rho)S + 2 r I_2 - (r + h + \xi + \sigma + \mu) I_1
```

for $`i>1`$:

``` math
\frac{\textstyle{dI_i}}{\textstyle{da}}  =  h(1-\rho)I_{i-1} + (i+1) r I_{i+1} - (i r + h + \xi + \sigma + \mu) I_1
```

## Implementation and Verification

In this implementation, the infinite system of equations is truncated at
$`N.`$ For the last differential equation $`I_N,`$ no new infections
occur ($`h=0`$) and there is no $`I_{N+1}`$ so nothing is added from
above:

``` math
\frac{\textstyle{dI_N}}{\textstyle{da}}  =  h(1-\rho)I_{N-1} - (N r + \xi + \sigma + \mu) I_N
```

In the implementation, if no value for $`N`$ is passed, it is set to

    N = round(max(10*h/r,20))

### Hybrid Variables:

This implementation also computes several other variables: Let
``` math
I = \sum_i I_i
```

The dynamics of the hybrid variable $`I`$ are:

``` math
\frac{\textstyle{dI}}{\textstyle{da}}  =  h(1-\rho)S - r I_1 - (\rho h+\xi+\sigma+\mu)I
```

Let $`m_n`$ denote the $`n^{th}`$ moment of the distribution of the MoI:

``` math
m_n = \frac{\textstyle{1}}{\textstyle{H}}\sum_i i^n  I_i
```

With some algebra (see below), we can derive an equation for the
dynamics of the hybrid variable $`m_1`$ are:
``` math
\frac{\textstyle{dm_1}}{\textstyle{da}}  =  h(1-\rho)\left(1-\frac{\textstyle{P}}{\textstyle{H}}\right) - \left(r + h \rho +\sigma + \xi\right) m_1
```

Similarly, we can compute the dynamics of $`m_2:`$

``` math
\frac{\textstyle{dm_2}}{\textstyle{da}}  =  h(1-\rho)\left( 1- \frac{\textstyle{P}}{\textstyle{H}} + 2 m_1\right)  - \left(2r + h \rho + \sigma + \xi \right) m_2 + r m_1
```

The dynamics of the hybrid variables, $`I,`$$`m_1,`$ and $`m_2`$ are
verified by solving the equations and using the formulas to compute
their values.

### Numerical Verification

``` r

library(ramp.falciparum)
library(ramp.func)
```

``` r

F_a = make_F_a(1)
q_out <- solve_SquIP(8/365, F_a, sigma=2/365, xi=0, Amax=3*365)
```

For verification, the variable $`m_1`$ is computed two ways, and we  
we note that they match. To see it, we plotted both, one in solid black
and the other dashed yellow:

``` r

with(q_out, plot(age, m_1, type = "l", ylab = expression(m[1]), main = "First Moment of the MoI"))
with(q_out, lines(age, m1, lty=2, col = "yellow"))
```

![](SquIP_files/figure-html/unnamed-chunk-3-1.png)

Since it’s hard to see differences, we can simply plot the differences.
Here, we’ve set the limits to be $`10^{-14}`$ so we can visualize the
differences:

``` r

ylm = 1e-14
with(q_out, plot(age, m_1-m1, type = "l", ylab = "Errors", main = "First Moment: Numerical Errors", ylim = c(-ylm, ylm)))
```

![](SquIP_files/figure-html/unnamed-chunk-4-1.png)

``` r

with(q_out, plot(age/365, m_2, type = "l", main = "Second Moment of the MoI",
                 ylab = expression(m[2]), 
                 xlab = "Age", ylim = range(m_2, m2)))
with(q_out, lines(age/365, m2, lty=2, col = "yellow"))
```

![](SquIP_files/figure-html/unnamed-chunk-5-1.png)

``` r

ylm = 1e-14
with(q_out, plot(age, m_2-m2, type = "l", 
                 ylab = "Errors", 
                 main = "Second Moment: Numerical Errors", 
                 ylim = c(-ylm, ylm)))
```

![](SquIP_files/figure-html/unnamed-chunk-6-1.png)

We note that if we had truncated the system of equations at $`N=4,`$ the
moments diverge:

``` r

q1_out <- solve_SquIP(8/365, F_a, sigma=2/365, xi=0, Amax=3*365, N=4)
with(q1_out, plot(age, m_1, type = "l", 
                  ylab = expression(m[1]), 
                  main = "First Moment of the MoI, N=4"))
with(q1_out, lines(age, m1, lty=2, col = "yellow"))
```

![](SquIP_files/figure-html/unnamed-chunk-7-1.png) \# Example

``` r

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
```

![](SquIP_files/figure-html/unnamed-chunk-8-1.png)

``` r

plot(tt, 0.1*F_t(tt), type = "l", lwd=2, ylim = c(0,0.04))
lines(aa, F_a(aa), type = "l", col = clrs[2])
lines(aa+365, F_a(aa, 365), type = "l", col = clrs[3])
lines(aa+730, F_a(aa, 730), type = "l", col = clrs[4])
lines(aa+1095, F_a(aa, 1095), type = "l", col = clrs[5])
lines(aa+1460, F_a(aa, 1460), type = "l", col = clrs[6])
```

![](SquIP_files/figure-html/unnamed-chunk-9-1.png)

``` r

plot(aa, F_a(aa), type = "n", lwd=2, ylim = c(0,0.04))
lines(aa, F_a(aa), type = "l", col = clrs[2])
lines(aa, F_a(aa, 365), type = "l", col = clrs[3])
lines(aa, F_a(aa, 730), type = "l", col = clrs[4])
lines(aa, F_a(aa, 1095), type = "l", col = clrs[5])
lines(aa, F_a(aa, 1460), type = "l", col = clrs[6])
```

![](SquIP_files/figure-html/unnamed-chunk-10-1.png)

``` r

plot(aa, cumsum(F_a(aa)), ylab = "Age", xlab = "Cumulative Exposure", type = "n", ylim = c(0,3))
lines(aa, cumsum(F_a(aa)), type = "l", col = clrs[2])
lines(aa, cumsum(F_a(aa, 365)), type = "l", col = clrs[3])
lines(aa, cumsum(F_a(aa, 730)), type = "l", col = clrs[4])
lines(aa, cumsum(F_a(aa, 1095)), type = "l", col = clrs[5])
lines(aa, cumsum(F_a(aa, 1460)), type = "l", col = clrs[6])
```

![](SquIP_files/figure-html/unnamed-chunk-11-1.png)

``` r

q3_out1 <- solve_SquIP(5/365, F_a, sigma=2/365, xi=0, Amax=5*365)
q3_out2 <- solve_SquIP(5/365, bday=365, F_a, sigma=2/365, xi=0, Amax=5*365)
q3_out3 <- solve_SquIP(5/365, bday=730, F_a, sigma=2/365, xi=0, Amax=5*365)
q3_out4 <- solve_SquIP(5/365, bday=365*3, F_a, sigma=2/365, xi=0, Amax=5*365)
q3_out5 <- solve_SquIP(5/365, bday=365*4, F_a, sigma=2/365, xi=0, Amax=5*365)
with(q3_out1, plot(age, m_1/x, col = clrs[2], ylab = "mean MoI", type ="l", ylim = c(1, 1.02)))
with(q3_out2, lines(age, m_1/x, col = clrs[3]))
with(q3_out3, lines(age, m_1/x, col = clrs[4]))
with(q3_out4, lines(age, m_1/x, col = clrs[5]))
with(q3_out5, lines(age, m_1/x, col = clrs[6]))
```

![](SquIP_files/figure-html/unnamed-chunk-12-1.png)

## Derivations

### Infected

We can derive the equation for $`dI/da`$ by collecting terms:

``` math
\frac{dI}{da} =\begin{array}{l} 
h(1-\rho)\left( S + \sum_i I_i \right) - h \sum_i I_i\\
\sum_i -ri I_i + (i+1) r I_{i+1} \\ 
- \sum_i \left(\sigma + \xi + \mu\right) I_i  \\
\end{array}
```

Since most things cancel, we simply get:

``` math
\frac{dI}{da} = h(1-\rho) S - r_1 I_1 - (\rho h + \xi + \sigma + \mu) I 
```

### MoI Distribution, First Moment

The dynamics of $`m_1`$ are given by:

``` math
\frac{\textstyle{dm_1}}{\textstyle{da}} = H^{-1} \sum_i i \frac{dI_i}{dt} - H^{-2} \sum_i i I_i \frac{d H}{da} 
```
Several of the terms can be handled all at once:
``` math
\frac 1 H
\sum_{i\geq 1}   i (\sigma + \xi + h + \mu) I_i = (\sigma + \xi + h + \mu) m_1  
```

**Natural Clearance** — For the terms involving $`r`$ – we can see a
pattern by looking at the terms one by one. For example, from
$`dI_2/dt`$ we add $`2 (3rI_3)`$, but from $`dI_3/dt,`$ we subtract
$`3(3rI_3)`$ leaving us with $`-3rI_3.`$ In general, we get:
``` math
\begin{array}{rl}
\frac{1}{H} \sum_{i\geq 1} i \left( -  i r I_i + (i+1) r I_{i+1} \right)
&= - \sum_i r i I_i/H \\
&= -  r m_1 
\end{array}
```

**Exposure** — For the terms involving exposure, $`h.`$ We note that
``` math
\begin{array}{l}
\frac{1}H \left[h (1-\rho) S  + \sum_{i\geq 1} (i+1) \left( h (1-\rho) I_i  \right) \right] \\ = h (1-\rho) \left(S + I + \sum_i i I_i \right) / H  \\ = h(1-\rho)\left(\frac{S+I}H + m_1\right) 
\end{array}
```

**Demography** —

Now, we compute $`dH/da`$ terms, we get:

``` math
- H^{-2} \sum_i iI_i \frac{dH}{da} = \frac{\sum_i i I_i \mu H}{H^2} = \mu m_1 
```
Putting it all together, we get that
``` math
\frac{\textstyle{dm_1}}{\textstyle{da}} = h(1-\rho) \frac{S+I}H - \left(r + \rho h + \sigma + \xi + \mu \right)  m_1  + \mu m_1
```
or
``` math
\begin{equation}
\frac{\textstyle{dm_1}}{\textstyle{da}} = 
h(1-\rho)\left(1 -  \frac{P}H \right) - \left(r + \rho h + \sigma + \xi  \right)  m_1  
\end{equation}
```

### MoI Distribution, Second Moment

The second, non-central moment of the distribution of the MoI ($`m_2`$)
is:
``` math
m_2 = \frac{\textstyle{\sum_i i^2 I_i}}H
```
As before, we compute:
``` math
\frac{\textstyle{dm_2}}{\textstyle{da}} = H^{-1} \sum_i i^2 \frac{dI_i}{dt} - H^{-2} \sum_i i^2 I_i \frac{d H}{da} 
```
We collect all the simple terms:
``` math
\frac{1}{H} \sum_{i\geq 1}   i^2 (\sigma + \xi + h + \mu) I_i = (\sigma + \xi + h + \mu) m_2
```

**Natural Clearance** — The terms involving $`r`$ are:
``` math
\begin{array}{l}
\frac{1}{H} \sum_{i\geq 1} i^2 \left( - i r I_i + (i+1) r I_{i+1} \right) \\
= \sum_i  r (i^3 - i(i-1)^2 I_i / H\\
= \sum_i  r (-2i^2 + 1) I_i / H\\
=  - 2r m_2 + r m 
\end{array}
```

**Exposure** — The terms with $`h (1-\rho)`$ are:

``` math
\begin{array}{c}
\frac{h (1-\rho)}H \left( S + \sum_{i\geq 1} (i+1)^2  I_i  \right) \\ 
= h(1-\rho) \left[ \left(\frac{S + I}H\right)  + 2m + m_2 \right] \\
\end{array}
```
and finally
``` math
- H^{-2} \sum_i i^2 I_i \frac{dH}{da} = \frac{\sum_i i^2 I_i \mu H}{H^2} = \mu m_2 
```
so we get
``` math
\frac{\textstyle{dm_2}}{\textstyle{da}} = h(1-\rho)\left(1-\frac{P}H + 2m\right) - \left(2r + \rho h + \sigma + \xi \right)  m_2 + rm 
```

# Random Variables for Parasite Infection Dynamics

We have derived mathematical formulas that describe the dynamics of
malaria infections as random variables in cohorts of humans as they age
(Henry JM, *et al.*, 2024)[^1] This R package is the computational
companion.

In the following, we review the mathematical formulas and the functions
in `ramp.falciparum` that compute the multiplicity of infection (MoI),
the age of infection (AoI), and the age of the youngest infection (AoY).

## Formulas and Computation

1.  Let $`z_d(\alpha,a)`$ denote the density for infections of age
    $`\alpha`$ in a host cohort of age $`a`$ born on day $`d`$.

    - The dynamics of $`z(\alpha,a)`$ are described by the following:
      ``` math
      \frac{\partial z}{\partial a} + \frac{\partial z}{\partial \alpha} = - r z
      ```
      with the boundary condition:
      ``` math
      \begin{equation}
       z_d(0,a) = h_d(a).
      \end{equation}
      ```
      Its solutions are given by:
      ``` math
      z_d(\alpha,a) = h_d(a-\alpha) e^{-r \alpha}
      ```

    - The function `zda` computes $`z_d(\alpha, a).`$

2.  The mean MoI is given by the formula:
    ``` math
    m_d(a) = \int_0^a z_d(\alpha, a) d \alpha
    ```

    - The distribution of the MoI is Poisson (Nåsell I, 1985)[^2] with
      mean $`m_d(a).`$

    - The function `meanMoI` computes $`m_d(a)`$ using `zda`

3.  The true prevalence is the
    ``` math
    p_d(a) = 1 - e^{-m_d(a)} 
    ```
    The function `truePR` computes the true prevalence using `meanMoI`

4.  The density function for the age of infection (AoI) is
    ``` math
    A_d(a) \sim f_A(\alpha, a, d) = \frac{z_d(\alpha,a)}{m_d(a)}
    ```
    and its moments are
    ``` math
    x_n(a, d | h) = \int_0^a  \alpha^n f_A(\alpha, a, d) d \alpha
    ```

5.  The age of the youngest infection (AoY) is defined as:
    ``` math
    Y_d(a) \sim f_Y(\alpha, a, d) = \min_{\zeta \sim M_d(a)} \left\{ \alpha_i \right\}_{i=1}^\zeta,  \alpha_i \sim A_d(a) 
    ```

    - The density function can be expressed in terms of the density and
      distribution functions of the AoI and MoI.
      ``` math
      f_Y(\alpha; a, d) = 
      f_A(\alpha, a,d) e^{-m_d(a) F_A(\alpha, a,d)} \frac{m_d(a)}{p_d(a)}.
      ```

    - The distribution function for the AoY is:
      ``` math
      F_Y(a) \sim
      \frac{1-e^{-m_d(a)F_A(\alpha, a,d)}}{1-e^{-m_d(a)}} = 
      \frac{1-e^{-m_d(a)F_A(\alpha, a,d)}}{p_d(a)} 
      \label{FY}
      ```

    - Its moments are:
      ``` math
      y_n(a, d | h) = \int_0^a  \alpha^n f_Y(\alpha | a, d, h) d \alpha 
      ```

6.  We also developed functions to compute the age of the youngest of
    $`N`$ infections, called YoN
    ``` math
    N_d(a) \sim
    \min_{N} \left\{ \alpha_i \right\}_{i=1}^N \mbox { where }  \alpha_i \sim A_d(a) 
    ```

    - The distribution function for YoN, $`N_d(a)`$, is
      ``` math
      F_N(\alpha, a, t) \sim 1- (1-F_A(\alpha, a, d))^N
      ```
      The following is a summary table of functions to compute the MoI,
      AoI, AoY, and all their moments.

    - The density function for YoN is found by differentiating:

``` math
f_N(\alpha, a, t) \sim N (1-F_A(\alpha, a, d))^{N-1}\frac{f_A(\alpha, a, d)}{m_d(a)}
```

## Quick Reference

The following is a summary table of functions to compute the MoI, AoI,
AoY, YoN, and all their moments.

|  | MoI | AoI | AoY | YoN |
|----|:--:|:--:|:--:|:--:|
|  | $`\zeta`$ | $`\alpha`$ | $`\alpha`$ | $`\alpha`$ |
|  | $`\zeta \geq 0`$ | $`0 \leq \alpha \leq a`$ | $`0 \leq \alpha \leq a`$ | $`0 \leq \alpha \leq a`$ |
| Random Variable | $`M_d(\zeta, a, h)`$ | $`A_d(\alpha, a, h)`$ | $`Y_d(\alpha, a , h)`$ | $`N_d(\alpha, a, h)`$ |
| Density Function | $`f_M(\zeta, a, h)`$ | $`f_A(\zeta, a, h)`$ | $`f_Y(\zeta, a, h)`$ | $`f_N(\zeta, a, h)`$ |
|  | dpois | dAoI | dAoY | dYoN |
| Distribution Function | $`F_M(\zeta, a, h)`$ | $`F_A(\zeta, a, h)`$ | $`F_Y(\zeta, a, h)`$ | $`F_N(\zeta, a, h)`$ |
|  | ppois | pAoI | pAoY | pYoN |
| Random Numbers | $`\hat M_d(\zeta, a, h)`$ | $`\hat A_d(\alpha, a, h)`$ | $`\hat Y_d(\alpha, a , h)`$ | $`\hat N_d(\alpha, a , h)`$ |
|  | rpois | rAoI | rAoY | rYoN |
| Moments | $`m_d(a, h)`$ | $`x_n(a, d, h)`$ | $`y_n(a, d, h)`$ |  |
|  | meanMoI | momentAoI | momentAoY |  |

## Demonstration

### Force of Infection (FoI)

``` r

clrs = viridisLite::turbo(7)
set.seed(234)
Sa = makepar_F_type2()
Sp = makepar_F_sin()

F_t <- make_ts_function(scale = 1, season_par=Sp)
FoI_a <- make_F_a(avg = 3/365, age_par=Sa, season_par=Sp)

tt <- seq(0, 3650, by=5)
aa <-seq(0, 365*5, by =5) 
plot(tt, 0.05*F_t(tt), type = "l")
```

![](Dynamics_files/figure-html/unnamed-chunk-2-1.png)

To compute anything, we must first set up a function to describe
exposure (see the
[FoI](https://dd-harp.github.io/ramp.falciparum/articles/FoI.md)
vignette). We define functions that plot the FoI for a cohort as it ages
(in red), but we can also compute the population average FoI (in black).
Different cohorts would experience different histories of exposure.

![](Dynamics_files/figure-html/unnamed-chunk-3-1.png)

### Computing `zda`

The function `ramp.falciparum::zda(alpha, a, FoIpar, ...)` uses the
formula in Eq. 1 to compute the density of parasite infections in a
cohort of humans as it ages.

Using `zda,` we can compute the density of parasites in a cohort of any
age without solving a full system of equations. Given a function
describing the FoI in the population, $`h(t)`$, and the cohort birthday,
$`d.`$

``` r

# devtools::load_all()
```

``` r

alpha = 60
a = 6*365
zda(60, 6*365, FoI_a) 
```

    ## [1] 0.008962569

The following computes the density of infections of every age in a
cohort of age 3.

``` r

zz = zda(a3years, max(a3years), FoI_a)
```

When we plot $`z_d(\alpha, a)`$, we note that as $`\alpha`$ grows
larger, the parasite cohort gets older. When we plot parasite cohorts by
age, time is going backwards on the x-axis.

![](Dynamics_files/figure-html/unnamed-chunk-6-1.png)

Now, we can imagine what `zda` would look like for several different
host cohorts at age three, but who were born at different times. In
effect, we are taking a snapshot of the cohorts at the same age, but at
different times.

The curves are different because the hosts were born at different
months, and they thus experienced different levels of exposure over the
first two years of life. Here the annual FoI is 5 infections, per
person, per year ($`\bar h = 5/365`$):

![](Dynamics_files/figure-html/unnamed-chunk-7-1.png)

## Multiplicity of Infection (MoI)

We define a random variable $`M`$ describing the multiplicity of
infection (MoI). The distribution of the MoI is Poisson (see the
[MoI](https://dd-harp.github.io/ramp.falciparum/articles/MoI.md)
vignette).

``` math
M_d(a) \sim f_M(\zeta; a, d) = \mbox{Pois}(m_d(a))
```

Since $`z_d(\alpha, a)`$ describes the density of *all* infections of
age $`\alpha`$ in a cohort of age $`a`$, the density of all infections
must be the MoI. Since $`0 \leq \alpha < a`$, it must be true that:

``` math
\begin{equation}
\tag{2}
m_d(a) = \int_0^a z_d(\alpha, a) d \alpha
\end{equation}
```

The function that computes $`m_d(a)`$ is called `meanMoI.`

``` r

mm = meanMoI(a3years, FoI_a, hhat=5/365)
```

Here, we plot the average MoI in the host cohort as it ages:

![](Dynamics_files/figure-html/unnamed-chunk-9-1.png)

## Age of Infection (AoI)

We define a random variable $`A_d(a)`$ that describes the age of
infection (AoI), which is given by the formula

``` math
A_d(a) \sim f_A(\alpha; a, d) = \frac{z_d(\alpha,a)}{m_d(a)}
```

### The Density Function, `dAoI`

We can compute $`A_d(a)`$ using the density function `dAoI`:

``` r

f_A = dAoI(a3years, max(a3years), FoI_a)
```

![](Dynamics_files/figure-html/unnamed-chunk-11-1.png)

Now, as we plot the distribution of the AoI in cohorts at age two, born
at different months (as we did above), we notice that the distributions
have changed shapes:

![](Dynamics_files/figure-html/unnamed-chunk-12-1.png)

### The Distribution Function, `pAoI`

The distribution function for $`A_d(a)`$ is:

``` math
F_A(a) \sim \int_0^\alpha f_A(\alpha; a, d) d\alpha
```

``` r

F_A = pAoI(a3years, max(a3years), FoI_a)
```

If our functions work correctly, then we should get approximately the
same answer from computing the cumulative sum of `dAoI.`

``` r

F_A_alt = cumsum(f_A)
```

We shouldn’t expect the answers to be exactly the same, but they should
be close, with the `pAoI` in black.

``` r

par(mar = c(5,4,1,1))

plot(a3years, F_A, type = "l", 
     xlab = "Parasite Cohort Age", 
     ylab = expression(1-F[X](alpha, a, bday)), lwd=3)

lines(a3years, F_A_alt, col = "red", lwd=2, lty =2)
```

![](Dynamics_files/figure-html/unnamed-chunk-15-1.png)

![](Dynamics_files/figure-html/unnamed-chunk-16-1.png)

### Random Numbers, `rAoI`

The function `rAoI` uses `pAoI` to generate random numbers from
$`F_A(\alpha)`$

``` r

rhx = rAoI(10000, 3*365, FoI_a)
```

A simple visual check computes the empirical CDF for the random variates
against $`F_A(\alpha)`$ computed using `pAoI`

``` r

par(mar = c(5,4,1,2))
plot(stats::ecdf(rhx), xlim = c(0,1095), cex=0.2, main = "", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(list(F[A](alpha), paste("ecdf"))))
lines(a3years, F_A, col = "red", lty = 2, lwd=2)
```

![](Dynamics_files/figure-html/unnamed-chunk-18-1.png)

We can also plot the distribution functions.

![](Dynamics_files/figure-html/unnamed-chunk-19-1.png)

### AoI Moments

Let $`x`$ denote the first moment of of $`A_d(a)`$:
``` math
x_d(a) = \left< A_d(a) \right> = \int_0^\infty \alpha \frac{z_d(\alpha, a)} {m_d(a)}
```

Similarly, we let $`x_d(a)[n]`$ denote the higher order moments of
$`A_d(a)`$:
``` math
x_{[n]}(a, d) = \int_0^\infty \alpha^n \frac{z_d(\alpha, a)} {m_d(a)}
```

``` r

moment1 = momentAoI(a3years, FoI_a)
moment2 = momentAoI(a3years, FoI_a, n=2)
moment3 = momentAoI(a3years, FoI_a, n=3)
```

The first three moments of the AoY plotted over time. In the top plot,
we’ve also plotted the $`n^{th}`$ root of the $`n^{th}`$ moment.

![](Dynamics_files/figure-html/unnamed-chunk-21-1.png)

## Age of the Youngest Infection (AoY)

We have derived a random variable $`Y_d(a)`$ describing the age of the
youngest infection (AoY). The density function for the AoY is:

``` math
Y_d(a) \sim f_Y(\alpha; a, d) = f_A(\alpha, a,d) e^{-m_d(a) F_A(\alpha, a,d)} \frac{m_d(a)}{p_d(a)}
```
The distribution function is:

``` math
F_Y(a) \sim
\frac{1-e^{-m_d(a)F_A(\alpha, a,d)}}{1-e^{-m_d(a)}} = 
\frac{1-e^{-m_d(a)F_A(\alpha, a,d)}}{p_d(a)} 
```

The derivations are found in a Suppplement to Henry JM, *et al.* (2024).

The mean AoY is:

``` math
\left< Y_d(a) \right> = \int_0^a \alpha \; f_Y(\alpha, a, d) \; d\alpha 
```

And the higher order moments for the AoY are:

``` math
\left< Y_d(a)^n \right> = \int_0^n \alpha^n \; f_y(\alpha, a, d) \; d\alpha 
```

### AoY Density, `dAoY`

The density function is computed with the function `dAoY.`

``` r

f_Y = dAoY(a3years, 3*365, FoI_a)
```

We can compare $`f_Y(\alpha)`$ (in black) to $`f_A(\alpha)`$ (in grey).

![](Dynamics_files/figure-html/unnamed-chunk-23-1.png)

### Random Variables, `rAoY`

``` r

raoy = rAoY(10^5, 3*365, FoI_a)
```

``` r

hist(raoy, breaks=seq(0, 1095, by = 15), 
     right=F, probability=T, main = "", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     border = grey(0.5)) -> out
lines(a3years, f_Y, type = "l", col = "red") 
```

![](Dynamics_files/figure-html/unnamed-chunk-25-1.png)

### AoY Moments

``` r

aa = seq(5, 3*365, by = 5) 
moment1y = momentAoY(aa, FoI_a)
moment2y = momentAoY(aa, FoI_a, n=2)
moment3y = momentAoY(aa, FoI_a, n=3)
```

The first three moments of the AoY plotted over time. In the top plot,
we’ve also plotted the $`n^{th}`$ root of the $`n^{th}`$ moment.

![](Dynamics_files/figure-html/unnamed-chunk-27-1.png)

------------------------------------------------------------------------

**Next:**

- In the vignette
  [MoI](https://dd-harp.github.io/ramp.falciparum/articles/MoI.md), we
  show that the mean MoI computed using the formula in Eq. 2 gives the
  same answer as other approaches.

[^1]: Henry JM, Carter AR, Wu SL, Smith DL (in preparation). A
    Probabilistic Synthesis of Malaria Epidemiology: Exposure,
    Infection, Parasite Densities, and Detection.

[^2]: Nåsell I (1985). Hybrid Models of Tropical Infections, 1st
    edition. Springer-Verlag.
    <https://doi.org/10.1007/978-3-662-01609-1>

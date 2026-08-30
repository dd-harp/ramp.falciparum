# Solve the queuing model \\M/M/\infty\\

A wrapper to solve the queuing model \\M/M/\infty\\ (see
[dMoIda](https://dd-harp.github.io/ramp.falciparum/reference/dMoIda.md)).

The function automatically sets the maximum MoI to be computed, and it
sets initial conditions. The equations are solved using
[deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html) and returned at
regular intervals dt from age 0 up to Amax (in days).

## Usage

``` r
solve_SIPm(
  h,
  FoI_a,
  bday = 0,
  r = 1/200,
  rho = 0.2,
  sigma = 1/365,
  xi = 1/365,
  eta = 1/25,
  mu = 0,
  H = 1000,
  Amax = 730,
  dt = 1
)
```

## Arguments

- h:

  the force of infection

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- r:

  the clearance rate for a simple infection

- rho:

  the fraction of incident cases that gets treted

- sigma:

  treatment rate for infected individuals

- xi:

  background drug taking

- eta:

  loss of chemoprotection

- mu:

  population death rate

- H:

  population size

- Amax:

  The maximum runtime (in days)

- dt:

  The output frequency (in days)

## Value

a [list](https://rdrr.io/r/base/list.html) with the orbits by name

## See also

[dMoIda](https://dd-harp.github.io/ramp.falciparum/reference/dMoIda.md)

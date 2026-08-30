# Solve the queuing model \\M/M/\infty\\

A wrapper to solve the queuing model \\M/M/\infty\\ (see
[dMoIda](https://dd-harp.github.io/ramp.falciparum/reference/dMoIda.md)).

The function automatically sets the maximum MoI to be computed, and it
sets initial conditions. The equations are solved using
[deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html) and returned at
regular intervals dt from age 0 up to Amax (in days).

## Usage

``` r
solveMMinfty(h, FoI_a, r = 1/200, bday = 0, Amax = 730, dt = 1)
```

## Arguments

- h:

  the force of infection

- FoI_a:

  a cohort trace function

- r:

  the clearance rate for a simple infection

- bday:

  the cohort birthday

- Amax:

  The maximum runtime (in days)

- dt:

  The output frequency (in days)

## Value

a [list](https://rdrr.io/r/base/list.html) with the orbits by name

## See also

[dMoIda](https://dd-harp.github.io/ramp.falciparum/reference/dMoIda.md)

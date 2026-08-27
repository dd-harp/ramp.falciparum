# Solve the system of differential equations to compute the approximate moments of the AoY over time.

Solve the system of differential equations to compute the approximate
moments of the AoY over time.

## Usage

``` r
solve_dAoYda(h, FoIpar, r = 1/200, tau = 0, Amax = 730, dt = 1, n = 8)
```

## Arguments

- h:

  the force of infection

- FoIpar:

  a FoI trace function

- r:

  the clearance rate for a simple infection

- tau:

  the cohort birthday

- Amax:

  The maximum runtime (in days)

- dt:

  The output frequency (in days)

- n:

  The number of terms to use in \\\phi(r,m)\\

## Value

a [data.frame](https://rdrr.io/r/base/data.frame.html) describing the
orbits

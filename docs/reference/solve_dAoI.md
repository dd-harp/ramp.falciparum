# Solve the system of differential equations to compute the moments of the AoI over time.

Solve the system of differential equations to compute the moments of the
AoI over time.

## Usage

``` r
solve_dAoI(h, FoIpar, r = 1/200, tau = 0, Amax = 730, dt = 1, N = 3)
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

- N:

  The total number of moments to compute

## Value

a data.frame with the orbits

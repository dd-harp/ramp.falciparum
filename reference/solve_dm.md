# Solve the hybrid model for the MoI

Solve the hybrid model for the MoI

## Usage

``` r
solve_dm(h, FoIpar, r = 1/200, tau = 0, Amax = 730, dt = 1)
```

## Arguments

- h:

  the average, annual force of infection

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

## Value

a [matrix](https://rdrr.io/r/base/matrix.html) with the orbits

# Solve the hybrid model for the MoI

Solve the hybrid model for the MoI

## Usage

``` r
solve_dm(h, FoI_a, r = 1/200, bday = 0, Amax = 730, dt = 1)
```

## Arguments

- h:

  the average, annual force of infection

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

a [matrix](https://rdrr.io/r/base/matrix.html) with the orbits

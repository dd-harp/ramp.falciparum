# Solve a system of differential equations to compute the true PR and the MoI

Solve a system of differential equations to compute the true PR and the
MoI

## Usage

``` r
solve_dpda(h, FoI_a, r = 1/200, bday = 0, Amax = 730, dt = 1)
```

## Arguments

- h:

  a scaling factor on the FoI

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

a [data.frame](https://rdrr.io/r/base/data.frame.html) with the orbits

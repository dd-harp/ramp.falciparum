# Random numbers for the AoI

Draw random numbers for the AoI from a cohort, \\\hat A\_\tau(a)\\

## Usage

``` r
rAoI(N, a, FoIpar, tau = 0, hhat = 1, r = 1/200, alphamin = 0)
```

## Arguments

- N:

  the number of observations

- a:

  host age

- FoIpar:

  a compound [list](https://rdrr.io/r/base/list.html) to compute
  \\h\_\tau(a)\\

- tau:

  the host cohort's birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

- alphamin:

  the minimum value of the AoI to return

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(alpha)

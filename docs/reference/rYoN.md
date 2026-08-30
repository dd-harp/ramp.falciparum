# The youngest of N infections, random numbers

The youngest of N infections, random numbers

## Usage

``` r
rYoN(R, N, a, FoI_a, bday = 0, hhat = NULL, r = 1/200, alphamin = 0)
```

## Arguments

- R:

  the number of observations

- N:

  the number of infections per person

- a:

  the host cohort age

- FoI_a:

  parameters that define an FoI function

- bday:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

- alphamin:

  the minimum value of the AoI to return

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length R

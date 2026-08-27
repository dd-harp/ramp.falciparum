# The youngest of N infections, distribution function

The youngest of N infections, distribution function

## Usage

``` r
pYoN(N, a, FoIpar, hhat = NULL, tau = 0, r = 1/200)
```

## Arguments

- N:

  the MoI

- a:

  the age of a cohort

- FoIpar:

  parameters that define an FoI function

- hhat:

  a local scaling parameter for the FoI

- tau:

  the cohort birthday

- r:

  the clearance rate of a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length a + 1

# The youngest of N infections, distribution function

The youngest of N infections, distribution function

## Usage

``` r
pYoN(N, a, FoI_a, bday = 0, hhat = NULL, r = 1/200)
```

## Arguments

- N:

  the MoI

- a:

  the age of a cohort

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate of a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length a + 1

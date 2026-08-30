# The youngest of N infections, density function

The youngest of N infections, density function

## Usage

``` r
dYoN(N, a, FoI_a, hhat = NULL, bday = 0, r = 1/200)
```

## Arguments

- N:

  the number of infections

- a:

  the age of a cohort

- FoI_a:

  a cohort trace function

- hhat:

  a local scaling parameter for the FoIs

- bday:

  the cohort birthday

- r:

  the clearance rate of a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length a + 1

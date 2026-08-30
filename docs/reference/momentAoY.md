# Compute the moments for the AoY density function for a cohort of age a

Compute the moments for the AoY density function for a cohort of age a

## Usage

``` r
momentAoY(a, FoI_a, bday = 0, hhat = 1, r = 1/200, n = 1)
```

## Arguments

- a:

  cohort age

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

- n:

  the moment desired

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

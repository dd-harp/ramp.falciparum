# Compute the moments for the AoY density function for a cohort of age a

Compute the moments for the AoY density function for a cohort of age a

## Usage

``` r
momentAoY(a, FoIpar, tau = 0, hhat = 1, r = 1/200, n = 1)
```

## Arguments

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

- n:

  the moment desired

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

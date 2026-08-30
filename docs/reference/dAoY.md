# The density function for the age of the youngest infection (AoY)

The density function for the age of the youngest infection (AoY)

## Usage

``` r
dAoY(alpha, a, FoI_a, bday = 0, hhat = 1, r = 1/200)
```

## Arguments

- alpha:

  the age of an infection, \\\alpha\\

- a:

  the age of the cohort

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- hhat:

  scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(alpha)

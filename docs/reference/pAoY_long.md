# Alternative method for computing the distribution function for the age of the youngest infection (AoY)

This method computes the distribution function for the AoY by summing
over AoI distributions

## Usage

``` r
pAoY_long(alpha, a, FoI_a, bday = 0, hhat = 1, r = 1/200)
```

## Arguments

- alpha:

  the age of an infection, \\\alpha\\

- a:

  cohort age

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

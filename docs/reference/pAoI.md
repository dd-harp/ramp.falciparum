# Compute the distribution function for AoI

The distribution function for the AoI is given by \$\$F_A(\alpha \| a,
\bday, h) = \int_0^\alpha f_A(\alpha, a, \bday \| h) d \alpha\$\$

## Usage

``` r
pAoI(alpha, a, FoI_a, bday = 0, hhat = 1, r = 1/200)
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

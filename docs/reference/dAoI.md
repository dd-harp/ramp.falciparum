# Compute the density function for AoI

The density of the AoI is given by \$\$f_A(\alpha \| a, \bday, h) =
\frac{\int_0^a z(\alpha, a)}{m\_\bday(a)}\$\$

## Usage

``` r
dAoI(alpha, a, FoI_a, bday = 0, hhat = 1, r = 1/200)
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

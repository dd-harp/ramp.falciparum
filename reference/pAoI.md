# Compute the distribution function for AoI

The distribution function for the AoI is given by \$\$F_A(\alpha \| a,
\tau, h) = \int_0^\alpha f_A(\alpha, a, \tau \| h) d \alpha\$\$

## Usage

``` r
pAoI(alpha, a, FoIpar, tau = 0, hhat = 1, r = 1/200)
```

## Arguments

- alpha:

  the age of an infection, \\\alpha\\

- a:

  host age

- FoIpar:

  a compound [list](https://rdrr.io/r/base/list.html) to compute
  \\h\_\tau(a)\\

- tau:

  the host cohort's birthday

- hhat:

  scaling parameter for
  [FoI](https://dd-harp.github.io/ramp.falciparum/reference/FoI.md)

- r:

  the clearance rate for a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(alpha)

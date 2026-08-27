# Compute the density function for AoI

The density of the AoI is given by \$\$f_A(\alpha \| a, \tau, h) =
\frac{\int_0^a z(\alpha, a)}{m\_\tau(a)}\$\$

## Usage

``` r
dAoI(alpha, a, FoIpar, tau = 0, hhat = 1, r = 1/200)
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

# Compute infection density in a cohort of humans, \\z\_\tau(\alpha, a \|h)\\

Given a function describing the FoI (\\h\_\tau(a)\\), and a parameter
describing the clearance rate of infections (\\r\\), the density of
parasite clones of age \\\alpha\\ distributed among a cohort of humans
of age \\a\\ is \$\$z\_\tau(\alpha, a) = e^{-r \alpha}
h\_\tau(a-\alpha)\$\$

## Usage

``` r
zda(alpha, a, FoIpar, tau = 0, hhat = 1, r = 1/200)
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

## See also

[FoI](https://dd-harp.github.io/ramp.falciparum/reference/FoI.md)

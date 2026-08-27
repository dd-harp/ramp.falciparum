# The mean MoI in a host cohort of age \\a\\

The mean multiplicity of infection (MoI) is \$\$m\_\tau(a\|h) = \int_0^a
z\_\tau(\alpha, a\|h) d \alpha\$\$

## Usage

``` r
meanMoI(a, FoIpar, tau = 0, hhat = 1, r = 1/200)
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

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

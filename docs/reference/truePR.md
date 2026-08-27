# Compute the true PR in a cohort as a function of age and exposure

The true PR is \$\$p\_\tau(a\|h) = 1 - e^{-m\_\tau(a\|h)}\$\$

## Usage

``` r
truePR(a, FoIpar, tau = 0, hhat = 1, r = 1/200)
```

## Arguments

- a:

  the host cohort age

- FoIpar:

  parameters that define an FoI function

- tau:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

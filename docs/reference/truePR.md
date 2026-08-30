# Compute the true PR in a cohort as a function of age and exposure

The true PR is \$\$p\_\bday(a\|h) = 1 - e^{-m\_\bday(a\|h)}\$\$

## Usage

``` r
truePR(a, FoI_a, bday = 0, hhat = 1, r = 1/200)
```

## Arguments

- a:

  the host cohort age

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

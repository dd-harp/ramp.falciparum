# Compute infection density in a cohort of humans, \\z\_\bday(\alpha, a \|h)\\

Given a function describing the FoI (\\h(a,d)\\), and a parameter
describing the clearance rate of infections (\\r\\), the density of
parasite clones of age \\\alpha\\ distributed among a cohort of humans
of age \\a\\ is \$\$z_d(\alpha, a) = e^{-r \alpha} h(a-\alpha, d)\$\$

## Usage

``` r
zda(alpha, a, FoI_a, bday = 0, hhat = 1, r = 1/200)
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

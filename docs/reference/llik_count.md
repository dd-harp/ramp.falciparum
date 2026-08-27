# Compute the likelihood for a single parasite count datum

Compute the likelihood for a single parasite count datum

## Usage

``` r
llik_count(count, mu, sz, par_O, q = 6, bvm = 13)
```

## Arguments

- count:

  the parasite count

- mu:

  the expected value for the count, given age

- sz:

  the size parameter for a negative binomial distribution

- par_O:

  parameters
  [d_Omega](https://dd-harp.github.io/ramp.falciparum/reference/d_Omega.md)

- q:

  the sampling volume (in \\\log\_{10} microliters\\)

- bvm:

  the blood volume (in \\\log\_{10} microliters\\)

## Value

the negative log likelihood for the observation

# Disribution function for the beta distribution, an alternative parameterization

Disribution function for the beta distribution, an alternative
parameterization

## Usage

``` r
pbeta1(p, mu, pSig = par_sigma_abc())
```

## Arguments

- p:

  a vector of probabilities

- mu:

  the mean value for the distribution (0 \<= mu \<= 1)

- pSig:

  parameters to dispatch the S3 function
  [sigma_mu](https://dd-harp.github.io/ramp.falciparum/reference/sigma_mu.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(p)

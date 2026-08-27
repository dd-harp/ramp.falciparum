# The random generation function for the beta distribution, an alternative parameterization Title

The random generation function for the beta distribution, an alternative
parameterization Title

## Usage

``` r
rbeta1(n, mu, pSig = par_sigma_abc())
```

## Arguments

- n:

  number of observations

- mu:

  the mean value for the distribution (0 \<= mu \<= 1)

- pSig:

  parameters to dispatch the S3 function
  [sigma_mu](https://dd-harp.github.io/ramp.falciparum/reference/sigma_mu.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length n

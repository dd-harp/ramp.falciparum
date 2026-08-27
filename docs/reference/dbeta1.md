# Density function for the beta distribution, an alternative parameterization

The beta distribution is parameterized using the mean and a function
`sigma_mu` that computes the variance as a function of the mean

## Usage

``` r
dbeta1(x, mu, pSig = par_sigma_abc())
```

## Arguments

- x:

  a vector of quantiles

- mu:

  the mean value for the distribution (0 \<= mu \<= 1)

- pSig:

  parameters to dispatch the S3 function
  [sigma_mu](https://dd-harp.github.io/ramp.falciparum/reference/sigma_mu.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(x)

# A function that returns constrained values of the variance for the beta distrution as a function of the mean.

A function that returns constrained values of the variance for the beta
distrution as a function of the mean.

## Usage

``` r
# S3 method for class 'abc'
sigma_mu(mu, par)
```

## Arguments

- mu:

  the mean value for the distribution (0 \<= mu \<= 1)

- par:

  parameters to dispatch and configure the instances

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(mu)

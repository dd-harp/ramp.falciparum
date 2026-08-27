# Modified beta distribution, density function

The beta distribution, parameterized by the mean and variance, modified
to return a number between 0 and bvm (the number of red blood cells)

## Usage

``` r
# S3 method for class 'beta'
d_Omega(xi, mu, bvm = 13, par_Omega = par_Omega_beta())
```

## Arguments

- xi:

  vector of quantiles for \\\log\_{10}\\ parasite densities

- mu:

  the expected value for \\\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\\log\_{10}\\ red blood cells

- par_Omega:

  parameters to compute parasite densities as a function of mu

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(xi)

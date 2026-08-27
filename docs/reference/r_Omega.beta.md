# Modified beta distribution, random numbers

The beta distribution, parameterized by the mean and variance, modified
to return a number between 0 and bvm (the number of red blood cells)

## Usage

``` r
# S3 method for class 'beta'
r_Omega(n, mu, bvm = 13, par_Omega = par_Omega_beta())
```

## Arguments

- n:

  the number of observations

- mu:

  the expected value for \\\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\\log\_{10}\\ red blood cells

- par_Omega:

  parameters defining parasite densities as a function of the mu

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(xi)

# The density function for parasite densities as a function of the mean

The density function for parasite densities as a function of the mean

## Usage

``` r
p_Omega(xi, mu, bvm = 13, par_Omega = par_Omega_beta())
```

## Arguments

- xi:

  a vector of probabilities for \\log\_{10}\\ parasite densities

- mu:

  the expected value for \\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\log\_{10}\\ red blood cells

- par_Omega:

  parameters defining parasite densities as a function of the mu

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(xi)

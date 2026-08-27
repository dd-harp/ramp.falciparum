# The density function for parasite densities in a simple malaria infection

The density function for parasite densities in a simple malaria
infection

## Usage

``` r
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

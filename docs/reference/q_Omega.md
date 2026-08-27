# The quantile function for parasite densities in a simple malaria infection

The quantile function for parasite densities in a simple malaria
infection

## Usage

``` r
q_Omega(xi, mu, bvm = 13, par_Omega = par_Omega_beta())
```

## Arguments

- xi:

  a vector of quantiles

- mu:

  the expected value for \\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\log\_{10}\\ red blood cells

- par_Omega:

  parameters defining parasite densities as a function of the mu

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(xi)

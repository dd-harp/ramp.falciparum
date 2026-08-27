# Random generation of parasite densities from a simple malaria infection

Random generation of parasite densities from a simple malaria infection

## Usage

``` r
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

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(n)

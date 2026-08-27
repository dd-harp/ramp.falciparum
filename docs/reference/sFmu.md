# Compute mean, expected parasite densities `mu` for multiple values of the age of infection `alpha`

Compute mean, expected parasite densities `mu` for multiple values of
the age of infection `alpha`

## Usage

``` r
sFmu(alpha, W, par_Fmu)
```

## Arguments

- alpha:

  a parasite's age of infection

- W:

  immune tracking variables

- par_Fmu:

  a [list](https://rdrr.io/r/base/list.html) that defines a model

## Value

mean log10 parasite densities, a
[numeric](https://rdrr.io/r/base/numeric.html) vector of length(alpha)

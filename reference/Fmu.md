# Compute mean, expected parasite densities `mu` as a function of the age of infection `alpha`

Compute mean, expected parasite densities `mu` as a function of the age
of infection `alpha`

## Usage

``` r
Fmu(alpha, W, par_Fmu)
```

## Arguments

- alpha:

  the age of a parasite infection

- W:

  the immune tracking variables

- par_Fmu:

  a [list](https://rdrr.io/r/base/list.html) that defines a model

## Value

mean log10 parasite densities, a
[numeric](https://rdrr.io/r/base/numeric.html) vector of length(alpha)

# Compute immune tracking variables as a function of host age and exposure

The function dispatches on `class(par)`

## Usage

``` r
Wda(a, FoI_a, bday = 0, hhat = 1, par = par_Wda_none())
```

## Arguments

- a:

  age of a host cohort

- FoI_a:

  a cohort trace function

- bday:

  cohort birthday

- hhat:

  overrides the value of hbar in par

- par:

  parameters in a [list](https://rdrr.io/r/base/list.html)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

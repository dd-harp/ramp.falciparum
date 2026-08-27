# Negative binomial PDF for raw parasite counts

Negative binomial PDF for raw parasite counts

## Usage

``` r
# S3 method for class 'nb'
d_counts(x, xi, bvm = 13, pars = par_nb())
```

## Arguments

- x:

  raw parasite counts

- xi:

  mean \\\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\\log\_{10}\\ red blood cells

- pars:

  parameters that define a detection function

## Value

binary detection result

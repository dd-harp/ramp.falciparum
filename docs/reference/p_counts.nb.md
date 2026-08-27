# Detection of infection given parasitemia

Detection of infection given parasitemia

## Usage

``` r
# S3 method for class 'nb'
p_counts(x, xi, bvm = 13, pars = par_nb())
```

## Arguments

- x:

  raw parasite counts

- xi:

  \\\log\_{10}\\ parasite densities

- bvm:

  blood volume as \\\log\_{10}\\ red blood cells

- pars:

  parameters that define a detection function

## Value

binary detection result

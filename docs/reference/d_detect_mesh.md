# Detection of infection given parasitemia

Detection of infection given parasitemia

## Usage

``` r
d_detect_mesh(xi_mesh, xi_density, bvm = 13, par_sample = par_nb())
```

## Arguments

- xi_mesh:

  a mesh over \\\log\_{10}\\ parasite densities

- xi_density:

  the density distribution for parasites over xi_mesh

- bvm:

  blood volume as \\\log\_{10}\\ red blood cells

- par_sample:

  parameters that define a detection function

## Value

the fraction of infected hosts that would test positive

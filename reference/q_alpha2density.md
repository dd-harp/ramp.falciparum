# The quantile function for parasite densities in a simple malaria infection of age alpha

The quantile function for parasite densities in a simple malaria
infection of age alpha

## Usage

``` r
q_alpha2density(
  x,
  alpha,
  W = 0,
  a = 0,
  par_RBC = par_lRBC_static(),
  par_Fmu = par_Fmu_base(),
  par_Omega = par_Omega_beta()
)
```

## Arguments

- x:

  a vector of quantiles

- alpha:

  the age of a parasite infection

- W:

  the immune tracking variables

- a:

  host cohort age

- par_RBC:

  parameters to compute
  [log10RBC](https://dd-harp.github.io/ramp.falciparum/reference/log10RBC.md)

- par_Fmu:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- par_Omega:

  parameters to compute
  [q_Omega](https://dd-harp.github.io/ramp.falciparum/reference/q_Omega.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(x)

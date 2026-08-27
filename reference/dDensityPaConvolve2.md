# The density function for the sum of two infections

The density function for the sum of two infections

## Usage

``` r
dDensityPaConvolve2(
  x,
  a,
  FoIpar,
  hhat = NULL,
  tau = 0,
  r = 1/200,
  par_RBC = par_lRBC_static(),
  par_Fmu = par_Fmu_base(),
  par_Omega = par_Omega_beta(),
  pWda = par_Wda_none()
)
```

## Arguments

- x:

  host cohort age

- a:

  host cohort age

- FoIpar:

  parameters that define an FoI function

- hhat:

  a local scaling parameter for the FoI

- tau:

  the cohort birthday

- r:

  the clearance rate for a simple infection

- par_RBC:

  parameters to compute
  [log10RBC](https://dd-harp.github.io/ramp.falciparum/reference/log10RBC.md)

- par_Fmu:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- par_Omega:

  parameters to compute
  [d_Omega](https://dd-harp.github.io/ramp.falciparum/reference/d_Omega.md)

- pWda:

  parameters to dispatch
  [Wda](https://dd-harp.github.io/ramp.falciparum/reference/Wda.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) value

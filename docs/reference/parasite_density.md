# Compute \\B\_\tau(a \| h)\\

This function first computes the PDF for simple infections,
[d_clone_density](https://dd-harp.github.io/ramp.falciparum/reference/d_clone_density.md),
and then computes the CDF for the full distribution, summed from MoI = 1
up to a large number using convolutions.

## Usage

``` r
parasite_density(
  meshX,
  a,
  FoIpar,
  tau = 0,
  hhat = 1,
  r = 1/200,
  RBC_par = par_lRBC_static(),
  Fmu_par = par_Fmu_base(),
  Omega_par = par_Omega_beta(),
  pWda = par_Wda_none()
)
```

## Arguments

- meshX:

  a mesh over parasite densities

- a:

  host cohort age

- FoIpar:

  parameters that define an FoI function

- tau:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

- RBC_par:

  parameters to compute
  [log10RBC](https://dd-harp.github.io/ramp.falciparum/reference/log10RBC.md)

- Fmu_par:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- Omega_par:

  parameters to compute
  [d_Omega](https://dd-harp.github.io/ramp.falciparum/reference/d_Omega.md)

- pWda:

  parameters to dispatch
  [Wda](https://dd-harp.github.io/ramp.falciparum/reference/Wda.md)

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length meshX

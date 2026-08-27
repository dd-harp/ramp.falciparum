# Random generation for M parasite densities in a host cohort with MoI

Random generation for M parasite densities in a host cohort with MoI

## Usage

``` r
rRda(
  M,
  R,
  a,
  FoIpar,
  tau = 0,
  hhat = 1,
  r = 1/200,
  alphamin = 7,
  RBC_par = par_lRBC_static(),
  Fmu_par = par_Fmu_base(),
  Omega_par = par_Omega_beta(),
  pWda = par_Wda_none()
)
```

## Arguments

- M:

  the MoI

- R:

  the number of observations

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

- alphamin:

  the minimum value of alpha allowed

- RBC_par:

  parameters to compute
  [log10RBC](https://dd-harp.github.io/ramp.falciparum/reference/log10RBC.md)

- Fmu_par:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- Omega_par:

  parameters defining parasite densities as a function of the mu

- pWda:

  parameters to dispatch
  [Wda](https://dd-harp.github.io/ramp.falciparum/reference/Wda.md)

## Value

a R by M [matrix](https://rdrr.io/r/base/matrix.html)

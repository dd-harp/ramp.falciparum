# Compute the mean parasite counts in simple infections

Compute the mean parasite counts in simple infections

## Usage

``` r
mean_clone_counts(
  a,
  FoI_a,
  dx = 0.1,
  hhat = 1,
  bday = 0,
  r = 1/200,
  Fmu_par = par_Fmu_base(),
  Omega_par = par_Omega_beta(),
  pWda = par_Wda_none(),
  par_sample = par_nb()
)
```

## Arguments

- a:

  host cohort age

- FoI_a:

  parameters that define an FoI function

- dx:

  the width of the mesh for computing the CDF of parasite densities

- hhat:

  a local scaling parameter for the FoI

- bday:

  the cohort birthday

- r:

  the clearance rate for a simple infection

- Fmu_par:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- Omega_par:

  parameters defining parasite densities as a function of the mu

- pWda:

  parameters to dispatch
  [Wda](https://dd-harp.github.io/ramp.falciparum/reference/Wda.md)

- par_sample:

  parameters that define a detection function

## Value

a [list](https://rdrr.io/r/base/list.html)

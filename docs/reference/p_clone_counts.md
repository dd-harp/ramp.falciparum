# Compute the distribution of parasite counts for simple infections

Compute the distribution of parasite counts for simple infections

## Usage

``` r
p_clone_counts(
  a,
  FoIpar,
  bins = c(1:5, 13),
  dx = 0.1,
  hhat = 1,
  tau = 0,
  r = 1/200,
  RBC_par = par_lRBC_static(),
  Fmu_par = par_Fmu_base(),
  Omega_par = par_Omega_beta(),
  pWda = par_Wda_none(),
  par_sample = par_nb()
)
```

## Arguments

- a:

  host cohort age

- FoIpar:

  parameters that define an FoI function

- bins:

  a set of break points for computing counts

- dx:

  the width of the mesh for computing the CDF of parasite densities

- hhat:

  a local scaling parameter for the FoI

- tau:

  the cohort birthday

- r:

  the clearance rate for a simple infection

- RBC_par:

  parameters to compute
  [log10RBC](https://dd-harp.github.io/ramp.falciparum/reference/log10RBC.md)

- Fmu_par:

  parameters to compute
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- Omega_par:

  parameters defining parasite densities as a function of mu

- pWda:

  parameters to dispatch
  [Wda](https://dd-harp.github.io/ramp.falciparum/reference/Wda.md)

- par_sample:

  parameters that define a detection function

## Value

a [list](https://rdrr.io/r/base/list.html)

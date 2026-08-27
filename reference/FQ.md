# The proportion of zero counts in

The proportion of zero counts in

## Usage

``` r
FQ(
  alpha,
  a,
  FoIpar,
  tau = 0,
  hhat = 1,
  r = 1/200,
  svm = 1e+06,
  par_RBC = par_lRBC_static(),
  par_Fmu = par_Fmu_base(),
  par_Omega = par_Omega_beta(),
  pWda = par_Wda_none(),
  par_sample = par_nb()
)
```

## Arguments

- alpha:

  the age of a parasite infection

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

- svm:

  the volume of a blood sample

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

- par_sample:

  parameters that define a detection function

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(xi)

# Compute \\P\_\tau(a \|h)\\

Compute \$\$P\_\tau(a\| h) \sim f_P(\xi; a, \tau \|h ) = \int_0^a
\Omega(\xi\|F\_\mu(\alpha)) \\ f_A(\alpha; a, \tau \| h) d\alpha\$\$

## Usage

``` r
d_clone_density(
  xi,
  a,
  FoIpar,
  tau = 0,
  hhat = 1,
  r = 1/200,
  par_RBC = par_lRBC_static(),
  par_Fmu = par_Fmu_base(),
  par_Omega = par_Omega_beta(),
  pWda = par_Wda_none()
)
```

## Arguments

- xi:

  log10

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

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(x)

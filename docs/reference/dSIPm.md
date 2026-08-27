# Derivatives for the SIP queuing model

This is the basic SIP queuing model, which starts from the model
\\M/M/\infty\\ but adds clearance through treatment and chemoprotection.
The model tracks the MoI in a cohort of humans as it ages. It assumes a
time- and age-dependent hazard rate for infection, called the force of
infection (FoI, \\h\_\tau(a)\\). Infections do not affect each other,
and each one clears independently at the rate \\r\\.

Let \\\zeta_i\\ the fraction of the population with MoI = i, then
\$\$\frac{d\zeta_0}{da}= -h\_\tau(a) \zeta_0 + r \zeta_1\$\$ and for
\\i\geq 1\\ \$\$\frac{d\zeta_i}{da}= h\_\tau(a) \left( \zeta\_{i-1} -
\zeta_i \right) - ri \zeta_i + r(i+1)\zeta\_{i+1}\$\$

This function computes the derivatives in a form that can be used by
[deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html).

## Usage

``` r
dSIPm(a, y, pars, FoIpar)
```

## Arguments

- a:

  the host age

- y:

  the state variables

- pars:

  the parameters

- FoIpar:

  \\h\_\tau(a)\\, a [list](https://rdrr.io/r/base/list.html) formatted
  to compute
  [FoI](https://dd-harp.github.io/ramp.falciparum/reference/FoI.md)

## Value

the derivatives as a [list](https://rdrr.io/r/base/list.html)

## See also

[solveMMinfty](https://dd-harp.github.io/ramp.falciparum/reference/solveMMinfty.md)

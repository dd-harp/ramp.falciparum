# Derivatives for [SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)

This computes the derivatives for the queuing model
[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md).
It also computes the hybrid variables \\I, m_1\\ and \\m_2\\

## Usage

``` r
dSquIPz(a, y, pars, FoIpar, Z)
```

## Arguments

- a:

  the host age

- y:

  a vector of state variables

- pars:

  the parameters

- FoIpar:

  \\h\_\tau(a)\\, a [list](https://rdrr.io/r/base/list.html) formatted
  to compute
  [FoI](https://dd-harp.github.io/ramp.falciparum/reference/FoI.md)

- Z:

  the MoE matrix

## Value

the derivatives as a [list](https://rdrr.io/r/base/list.html)

## See also

[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md) &
[solve_SquIP](https://dd-harp.github.io/ramp.falciparum/reference/solve_SquIP.md)

This function computes the derivatives in a form that can be used by
[deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html).

[solveMMinfty](https://dd-harp.github.io/ramp.falciparum/reference/solveMMinfty.md)

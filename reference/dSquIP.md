# Derivatives for [SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)

This function computes the derivatives for the queuing model
[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md),
in a form that can be used by
[deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html).

It also computes the derivatives for the hybrid variables \\m_1\\ and
\\m_2\\

## Usage

``` r
dSquIP(a, y, pars, FoIpar)
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

## Value

the derivatives as a [list](https://rdrr.io/r/base/list.html)

## See also

[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md) &
[solve_SquIP](https://dd-harp.github.io/ramp.falciparum/reference/solve_SquIP.md)

[solveMMinfty](https://dd-harp.github.io/ramp.falciparum/reference/solveMMinfty.md)

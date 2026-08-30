# Solve [SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)

This function solves the model
[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)
using [deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html). It is a
wrapper around the derivatives function
[dSquIP](https://dd-harp.github.io/ramp.falciparum/reference/dSquIP.md).
It does the following:

- `N` — set a value for the truncation parameter, if no value was passed

- `y` — set up the the initial values vector as a named vector

- `times` — set up a mesh over the independent variable \\a\\ (age, in
  days)

- `parms` — set up a named vector with the parameter values

After solving the equations, it parses the outputs, and returns the
values of the dependent variables by name. The wrapper

## Usage

``` r
solve_SquIP(
  h,
  F_a,
  bday = 0,
  r = 1/200,
  rho = 0.2,
  sigma = 1/365,
  xi = 1/365,
  eta = 1/25,
  mu = 0,
  H = 1000,
  Amax = 730,
  da = 1,
  N = NULL
)
```

## Arguments

- h:

  the force of infection

- F_a:

  a trace function

- bday:

  the cohort birthday

- r:

  the clearance rate for a simple infection

- rho:

  the fraction of incident cases that gets treted

- sigma:

  treatment rate for infected individuals

- xi:

  background drug taking

- eta:

  loss of chemoprotection

- mu:

  population death rate

- H:

  cohort size

- Amax:

  The maximum runtime (in days)

- da:

  The output interval (in days)

- N:

  truncation parameter (maximum MoI): if `NULL`, then set by rule

## Value

a named [list](https://rdrr.io/r/base/list.html) of parsed outputs

## Parased Outputs

After solving the equations, the solutions are parsed and returned as a
named list:

- `age` — the ages at which the dependent variables were output

- `S` — Susceptible, a vector

- `P` — Chemoprotected, a vector

- `Ii` — Infected \\\times\\ MoI, a matrix

- `H` — Cohort population size

- `m_1` — The hybrid variable \\m_1\\, mean MoI

- `m_2` — The hybrid variable \\m_2\\, the second moment of the MoI

- `m1`— The mean MoI, computed from `Ii`

- `m2`— The second moment of the MoI, computed from `Ii`

- `x` — The prevalence of infection: \\x = (H-P-S)/H\\

- `F_r` — The rate of natural clearance: \\r I_1\\

- `N` — The truncation parameter

- `out` — The matrix returned by `deSolve`

## See also

[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)

# SquIP

A model for the distribution of the multiplicity of infection (MoI) in a
cohort of humans as it ages. The model extends a basic queuing model
that has been used in malaria epidemiology before called \\M/M/\infty\\:
it adds parameters to model treatment with anti-malarial drugs, which
cures infections and is followed by a short period of chemo-protection.
The model is named after its states:

## **SquIP** \\=\\

**S** - **S**usceptible \\+\\

**quI** - a **qu**euing process for **I**nfected \\+\\

**P** - chemo-**P**rotected

The master equations are presented below. The software implementation of
**SquIP** has two parts:

- [dSquIP](https://dd-harp.github.io/ramp.falciparum/reference/dSquIP.md)
  computes the derivatives

- [solve_SquIP](https://dd-harp.github.io/ramp.falciparum/reference/solve_SquIP.md)
  is a wrapper around
  [dSquIP](https://dd-harp.github.io/ramp.falciparum/reference/dSquIP.md)
  that

  - sets the value of a truncation parameter

  - sets up the initial values, times, and parameters

  - calls [deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html) to
    solve the system of ordinary differential equations

  - parses the outputs

This software implementation truncates the infinite system of
dIifferential equations. To check the accuracy of the truncated system,
we have derived hybrid variables.

## The Dynamical System

**Variables**

The independent variable is age, \\a.\\

The dependent variables are:

- \\I_i\\ is the expected number infected with \\i\\ clones

- \\S\\ is the expected number uninfected and susceptible to infection
  \\(S=I_0)\\

- \\P\\ is the expected number uninfected and chemo-protected

- \\H = S+P+\sum_i I_i\\ is expected number in the cohort

The model assumes everyone is born susceptible so \\H\\ is passed as a
parameter, and the initial conditions are set to \\S(0) = H\\ and
\\I_i(0) = P(0)=0.\\

**Exposure** — \\h(a)\\

Let \\h(a)\\ denote the force of infection (FoI). This model assumes
that each incident infection increases the MoI by exactly one. The FoI
is passed as a trace function.

**Natural Clearance** — \\r\\

The model assumes that each clone clears independently of the others at
rate \\r,\\ so if a person has \\i\\ clones, the MoI goes down by one at
the rate \\ri.\\

**Treatment and Chemoprotection** — \\\xi, \rho, \sigma, \eta\\

These equations assume that individuals are treated and cured for
several different reasons:

- Everyone in the population who is not chemoprotected takes drugs at a
  background rate \\\xi\\

- Some incident infections cause disease, and a fraction \\\rho\\ gets
  treated

- In addition to these other modes of treatment, infected individuals
  get treated at the higher rate \\\sigma\\

After being treated and cured, individuals enter the chemoprotected
class \\P.\\ Chemoprotection is lost at the rate \\\eta,\\ and they
enter the susceptible class \\S.\\

**Demography** — \\\mu\\

The model assumes that all individuals die at the rate \\\mu(a)\\, such
that

\$\$\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H\$\$

## Master Equations

The dynamics are described by an infinite set of differential equations:
\$\$ \begin{array}{rl} dS/da &= \eta P + r I_1 - (h + \xi + \mu) S \\
dP/da &= \left(h \rho + \xi \right) \left(H-P \right) + \sigma \sum_i
I_i - \left(\eta + \mu \right) P \\ dI_1/ da &= h(1-\rho)S + 2 r I_2 -
(r + h + \xi + \sigma + \mu) I_1 \\ dI_i/da &= h(1-\rho)I\_{i-1} + (i+1)
r I\_{i+1} - (i r + h + \xi + \sigma + \mu) I_i\\ \end{array} \$\$

In this software implementation, the infinite system of equations is
truncated after \\i=N\\ and \$\$dI_N/da = h(1-\rho)I\_{N-1} - (N r + h +
\xi + \sigma + \mu) I_N\$\$

## Hybrid Variables

To check the accuracy of the truncated system, we have derived hybrid
variables for the first two moments of the distribution of the MoI.

Let \\m_n\\ denote the \\n^{th}\\ moment of the MoI: \$\$m_n =
\frac{\textstyle{1}}{\textstyle{H}}\sum_i i^n I_i\$\$

This software implementation computes the dynamics of the hybrid
variables \\m_1\\ and \\m_2\\: \$\$ \begin{array}{rl} d m_1 /da &=
h(1-\rho)\left(1-P/H\right) - \left(r + h \rho +\sigma + \xi\right) m_1
\\ d m_2 /da &= h(1-\rho)\left( 1- P/H + 2 m_1\right) - \left(2r + h
\rho + \sigma + \xi \right) m_2 + r m_1 \\ \end{array} \$\$

After solving, the output is parsed and the first two moments of the
distribution of the MoI are computed from \\I_i\\. The empirical moments
should match \\m_1\\ and \\m_2\\ unless the value of the truncation
parameter \\N\\ is set too low.

## See also

[dSquIP](https://dd-harp.github.io/ramp.falciparum/reference/dSquIP.md)
&
[solve_SquIP](https://dd-harp.github.io/ramp.falciparum/reference/solve_SquIP.md)
\| The model
[SquIPz](https://dd-harp.github.io/ramp.falciparum/reference/SquIPz.md)
extends **SquIP.** \|
[SIPm](https://dd-harp.github.io/ramp.falciparum/reference/SIPm.md) is
an approximating model.

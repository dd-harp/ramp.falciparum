# SquIPz

This is a model for the dIistribution of the multiplicity of infection
(MoI) in a cohort of humans as it ages. It extends the queuing model
[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md) by
adding a Tweedie process to model heterogeneous exposure: each infection
event increases the MoI by the multiplicity of exposure (MoE), which is
given by a probability mass function \\F\_\mbox{MoE}.\\ The model is
named after its states and processes:

## **SquIPz** \\=\\

**S** - **S**usceptible \\+\\

**quI** - a **qu**euing process for **I**nfected \\+\\

**P** - chemo-**P**rotected

**z** - the Tweedie process distribution

This software implementation of **SquIP** has two parts:

- [dSquIPz](https://dd-harp.github.io/ramp.falciparum/reference/dSquIPz.md)
  computes the derivatives

- [solve_SquIPz](https://dd-harp.github.io/ramp.falciparum/reference/solve_SquIPz.md)
  is a wrapper around
  [dSquIPz](https://dd-harp.github.io/ramp.falciparum/reference/dSquIPz.md)
  that

  - sets the value of a truncation parameter

  - sets up the initial values, times, and parameters

  - calls [deSolve::ode](https://rdrr.io/pkg/deSolve/man/ode.html) to
    solve the system of ordinary differential equations

  - parses the outputs

This software implementation truncates the infinite system of
dIifferential equations. To check the accuracy of the truncated system,
we have derived hybrid variables.

## The Model

**Variables** — The independent variable is age, \\a.\\

The dependent variables are:

- \\I_i\\ is the expected number infected with \\i\\ clones

- \\S\\ is the expected number uninfected and susceptible to infection
  \\(S=I_0)\\

- \\P\\ is the expected number uninfected and chemo-protected

- \\H = S+P+\sum_i I_i\\ is expected number in the cohort

The model assumes everyone is born susceptible so \\H\\ is passed as a
parameter, and the initial conditions are set to \\S(0) = H\\ and
\\I_i(0) = P(0)=0.\\

**Exposure** (\\h\\) —

Let \\h\\ denote the force of infection (FoI). This model assumes that
each incident infection increases the MoI by \\i\\ with probability
\$\$F\_\mbox{MoE}(X=i) = z_i.\$\$

**Natural Clearance** — The model assumes that each clone clears
independently of the others at rate \\r,\\ so if a person has \\i\\
clones, the MoI goes down by one at the rate \\ri.\\

**Treatment and Chemoprotection** — These equations assume that
individuals are treated and cured for several different reasons:

- Everyone in the population takes drugs at a background rate \\\xi\\

- Incident infections cause disease and a fraction \\\rho\\ gets treated

- In addition to other modes of treatment, infected individuals get
  treated at the higher rate \\\sigma\\

After being treated and cured, individuals enter the chemoprotected
class \\P.\\ Chemoprotection is lost at the rate \\\eta,\\ and they
enter the susceptible class \\S.\\

**Demography** – The model assumes that all individuals die at the rate
\\\mu(a)\\, such that

\$\$\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H\$\$

**Dynamics** — The dynamics are described by an infinite set of
differential equations: \$\$ \begin{array}{rl} dS/da &= \eta P + r I_1 -
(h + \xi + \mu) S \\ dP/da &= \left(h \rho + \xi \right) \left(H-P
\right) + \sigma \sum_i I_i - \left(\eta + \mu \right) P \\ dI_i / da
&=h (1-\rho) \left( S z_i + \sum_j I\_{i-j} z_j \right) + (i+1) r
I\_{i+1} - \left(h + \sigma + r i + \xi + \mu \right) I_i\\ \end{array}
\$\$

In this software implementation, the infinite system of equations is
truncated at \\N\\ and \$\$dI_N/da = h (1-\rho) \sum\_{i \geq N} \left(
S z_i + \sum_j I\_{i-j} z_j \right) - \left(h + \sigma + r N + \xi + \mu
\right) I_N\$\$

## Hybrid Variables

Let \\m_n\\ denote the \\n^{th}\\ moment of the MoI: \$\$m_n =
\frac{\textstyle{1}}{\textstyle{H}}\sum_i i^n I_i\$\$ and let \\\tilde
z_n\\ denote the moments of \\F\_\mbox{MoE}:\\ \$\$\tilde z_n = \sum_i
i^n z_i\$\$

This software implementation computes the dynamics of the hybrid
variables \\m_1\\ and \\m_2\\: \$\$ \begin{array}{rl} d m_1 /da &=
h(1-\rho)\tilde z_1 \left(1-P/H \right) - \left(r + h \rho +\sigma +
\xi\right) m_1 \\ d m_2 /da &= h(1-\rho)\hat z_2 \left( 1-P/H\right) +
(r+2 h(1-\rho)\hat z) m - \left(2r + h \rho + \sigma + \xi \right) m_2
\\ \end{array} \$\$ After solving, the output is parsed and the first
two moments of the distribution of the MoI are computed from \\I_i\\.
The empirical moments should match \\m_1\\ and \\m_2\\ unless the value
of the truncation parameter \\N\\ is set too low.

## See also

[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md)
and [SIPm](https://dd-harp.github.io/ramp.falciparum/reference/SIPm.md)

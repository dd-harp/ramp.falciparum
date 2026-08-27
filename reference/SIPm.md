# SIPm

**SIPm** is an enhanced compartmental model that approximates the model
**SquIPz**

## Variables

The independent variable is age, \\a.\\

Three variables make up a compartmental model:

- \\S\\ is the expected number uninfected and susceptible to infection

- \\P\\ is the expected number uninfected and chemo-protected

- \\I\\ is the expected number infected

Two variables track the distribution of the MoI

- \\m_1\\ is the first moment of the MoI

- \\m_2\\ is the second moment of the MoI

## Demography

We let \$\$H = S+I+P\$\$ denote the expected population density of the
cohort as it ages. The model assumes that all individuals die at the
rate \\\mu\\, such that

\$\$\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H\$\$

In the model, we assume that everyone is born susceptible so \\H\\ is
passed as a parameter, and the initial conditions are set to \\S(0) =
H\\ and \\I_i(0) = P(0)=0.\\

## Exposure

We let \\h\\ denote the force of infection (FoI), and we let
\\F\_\mbox{MoE}\\ denote the distribution of the multiplicity of
exposure (MoE).

- \\\tilde z_1\\ is the first moment of \\F\_\mbox{MoE}\\

- \\\tilde z_2\\ is the second moment of \\F\_\mbox{MoE}\\

## Clearance

In the queuing model
[SquIP](https://dd-harp.github.io/ramp.falciparum/reference/SquIP.md),
infections clear at the rate \\r I_1.\\ If we think of **SIPm** as an
approximating model, then \$\$r I_1 \approx F_r(m_1, m_2) I\$\$

To get \$F_r,\$ we assume that the distribution of the MoI in these
models is approximately a zero-inflated negatively binomial distributed.
We assume \$\$ \frac{I_i} H\sim \mbox{NB}\left(\zeta=i\| m, \phi
\right)\left(1-\frac PH\right) \$\$ where
`mu1=\eqn{m_1} is the mean and `size\`\\=\phi\\ (the size parameter) is
given by: \$\$\phi = \frac{m_1^2}{m_2-m_1^2-m_1}\$\$ so \$\$ F_r(m,
\phi) = \begin{cases} r & \mbox{if } m = 0 \\ r \frac{
\mbox{NB}\left(\zeta=1\| m, \phi \right) }{1- \mbox{NB}\left(\zeta=0 \|
m, \phi \right) } \left(1 - \frac{P}{H} \right) & \mbox{if } m\>0
\end{cases} \$\$

In effect, the approximation assumes \\I_1/H\\ is close to the
probability the MoI is equal to 1 in a zero-inflated negative binomial
distribution. The distributions of the MoI are probably not exactly
negatively binomially distributed. We use the Kullback-Liebler
divergence to evaluate how well a negative binomial approximates the
distributions in the queueing models.

## Treatment and Chemprotection

These equations assume that individuals are treated and cured for
several different reasons:

- Incident infections cause disease and a fraction \\\rho\\ gets treated

- Everyone in the population takes drugs at a background rate \\\xi\\

- In addition to other modes of treatment, infected individuals get
  treated at the higher rate \\\sigma\\

After being treated and cured, individuals enter the chemoprotected
class \\P.\\ Chemoprotection is lost at the rate \\\eta,\\ and they
enter the susceptible class \\S.\\

## Differential Equations

**Hybrid Variables** The compartmental dynamics are enhanced by
computing a hybrid variable, a non-compartmental state variable
describing the mean MoI (\\m_1\\), with dynamics: \$\$
\frac{\textstyle{dm_1}}{\textstyle{da}} = h (1-\rho) \hat z \left(
1-\frac{\textstyle{P}}{\textstyle{H}} \right) - \left(r + h \rho +
\sigma + \xi \right) m_1\$\$

We also compute the dynamics of the second moment of the distribution of
the MoI: \$\$ \begin{array}{rl} \frac{\textstyle{dm_2}}{\textstyle{da}}
= & h(1-\rho)\hat z_2 \left( 1- \frac{
\textstyle{P}}{\textstyle{H}}\right) + (r+2 h(1-\rho)\hat z) m_1
\\\[6pt\] &- \left(2r + h \rho + \sigma + \xi \right) m_2 \\
\end{array}\$\$

**Compartmental Variables**

The infinite infection states in the family of queuing models specified
by
[SquIPz](https://dd-harp.github.io/ramp.falciparum/reference/SquIPz.md)
are compressed down to a finite set of mutually exclusive and
collectively exhaustive states – \\S\\, \\I\\, and \\P\\ – with
dynamics: \$\$ \begin{array}{rl} \frac{\textstyle{dS}}{\textstyle{da}}
&= F_r(m, \phi) I + \eta P - (h + \xi + \mu) S \\\[6pt\]
\frac{\textstyle{dI}}{\textstyle{da}} &= h \left(1-\rho\right) S
\\\[6pt\] & - \left( \rho h + \xi + \sigma + F_r(m, \phi) + \mu \right)
I \\\[6pt\] \frac{\textstyle{dP}}{\textstyle{da}} &= (h \rho+
\xi)(H-P) + \sigma I - (\eta+\mu) P \\\[6pt\] \end{array}\$\$

## Parameters

## See also

[SquIPz](https://dd-harp.github.io/ramp.falciparum/reference/SquIPz.md)
and SIPm

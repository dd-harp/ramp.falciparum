# `ramp.falciparum`

## Falciparum Malaria Mathematical Epidemiology

## Install

To install the latest version of `ramp.falciparum` from GitHub, run the
following lines of code in an R session. We also recommend loading the
supporting package `ramp.func`

    library(devtools)
    devtools::install_github("dd-harp/ramp.falciparum")
    devtools::install_github("dd-harp/ramp.func")

The packages can then be loaded in the normal way:

    library(ramp.falciparum)
    library(ramp.func)

## Overview

The burden of falciparum malaria – human disease caused by infection
with *Plasmodium falciparum* – on human health and health systems is
enormous. In countries where malaria is endemic, the prevalence of
malaria varies, reflecting differences in exposure. In some places, most
people test positive and more than half of all outpatients get diagnosed
with malaria. Severe malaria is uncommon, but it can result in death.
Mild disease – subjective fever or malaise – is quite common.

There is a need to understand malaria epidemiology well enough to manage
it, but the epidemiology of falciparum malaria is complex, and a
synthesis has proven elusive. The management of malaria has two main
goals: there is a need to implement policies that **reduce the burden**
of disease; and most countries have the long term goal of **eliminating
malaria.** To manage malaria, we must measure it. Mathematical models
are one useful way of dealing with the complexity.

**`ramp.falciparum`** is a software project focused on the mathematical
epidemiology of falciparum malaria. The software has taken a
computational approach — the models have been encoded, and wrapper
functions were written to make it easy to solve them using
[**`deSolve`**](https://cran.r-project.org/web/packages/deSolve/index.html).
Each model has been described in a vignette, which also has any
mathematical derivations (*e.g.* for hybrid variables).

The software was written to make the mathematical epidemiology of
malaria more accessible to everyone, and to explain some recent
innovations. In describing the complexity of malaria without resorting
to individual-based simulation, we use **random variables.** In forging
a synthesis, we make extensive use of **hybrid variables** describing
the distributions of the age of infection (AoI), the multiplicity of
infection (MoI), and cumulative exposure. The hybrid variables are used
in **enhanced compartmental models:** a traditional compartmental model
is enhanced with hybrid variables.

## Malaria Epidemiology as Ontogeny

In developing this approach to malaria, we consider the developing of
immunity to malaria in a cohort as it ages over a lifetime The
mathematical formulations can be reduced to a simple question:

> **Given a history of exposure to malaria in a person of some age, what
> are the factors that determine the outcomes of the next infection?**

With this narrow focus, all the models track [**cohort
dynamics:**](https://dd-harp.github.io/ramp.falciparum/articles/Cohorts.md).
the models examine malaria in a cohort of humans as it ages. the
independent variable is thus *age* and not time. The force of infection
in these models (denoted $`h`$) is a [**trace
function**](https://dd-harp.github.io/ramp.falciparum/articles/Exposure.md).

Generically, let $`\mathbf{X}`$ denote a state space. Dynamical systems
describing malaria dynamics in a cohort are forced by a trace function
$`h(a,d)`$, the force of infection for a cohort at age $`a`$ born on day
$`d`$. The models in **`ramp.falciparum`** take the form:

``` math
\frac{\textstyle{d \mathbf{X(a)}}}{\textstyle{da}} =  F_{\textbf{X}}\left(h, \textbf{X}\right)
```

In malaria analytics, the focus is on the outcomes in populations.
Models for malaria analytics are found in a suite of related software
packages, collectively known as
[**SimBA.**](https://faculty.washington.edu/smitdave/simba/) In SimBA,
models track the temporal dynamics of malaria with age-structure (see
the [Cohort
Dynamics](https://dd-harp.github.io/ramp.falciparum/articles/Cohorts.md)
vignette).

## Related Software

**`ramp.falciparum`** is focused malaria epidemiology defined narrowly
to include exposure, infection, infectiousness, disease, immunity, and
diagnostics and detection.

It does not include malaria transmission dynamics and control or
mosquito ecology. Those topics are covered in related projects:

- Software developed to support malaria analytics, including malaria
  transmission dynamics and control, has been developed around
  [**`ramp.xds`**](https://dd-harp.github.io/ramp.xds/) and its
  satellite packages, collectively called
  [**SimBA**](https://faculty.washington.edu/smitdave/simba/). Of these,
  two are worth mentioning:

  - **`ramp.library`** holds a library of malaria models that can be
    used by **`ramp.xds`**

  - **`ramp.demog`** has utilities to develop age-structured models.

- A website devoted to [**Malaria
  Theory**](https://faculty.washington.edu/smitdave/malaria_theory/)
  makes extensive use of this software.

## Towards a Synthesis

To be useful in policy the models must be fit for purpose, realistic
enough to be compelling, accurate enough to useful, and simple enough to
be insightful. The software to support malaria analytics should support
**nimble model building:** it should be easy to build models, modify the
models as policy discussions evolve, and then run most analyses without
a supercomputing cluster.

Published models of malaria tend to focus on isolated aspects of
malaria. Models that include all the required features tend to be
within-host models, which are neither nimble nor simple enough to be
insightful. Using enhanced compartmental models, we can achieve a new
synthesis.

If we want to use models effectively in malaria policy, we will need to
go one step further and show that are models are realistic and accurate
enough for the task at hand. These issues fall under the category of
[**robust analytics for malaria
policy.**](https://dd-harp.github.io/ramp.falciparum/articles/RobustAnalytics.md)

## References

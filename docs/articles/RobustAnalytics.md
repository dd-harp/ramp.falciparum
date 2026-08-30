# Robust Analytics

Policy advice should be robust to uncertainty. To build rigor around the
idea of robustness in policy, we recognized the need for a bespoke
inferential framework designed for malaria analytics, which gave rise to
**RAMP** (Robust Analytics for Malaria Policy). An important goal of
RAMP was to build an integrated computational environment to support a
range of activities, including conventional statistical analysis and
methods to characterize, quantify, and propagate uncertainty. As RAMP
took shape, RAMP software was being developed to implemented the
principles in a stable computational form to facilitate the
transformation of data into robust policy advice. Today, core RAMP
software includes several R packages (see
[ramp.malaria](https://dd-harp.github.io/ramp.malaria/)) to deal with
malaria epidemiology, transmission dynamics, and control in the broad
sense. The [ramp.library](https://dd-harp.github.io/ramp.library/)
includes reusable code for one of several models of human malaria
infection and immuno-epidemiology.

We recognized that the complex epidemiology of malaria - defined in the
narrow sense to include only parasite infection in humans and processes
affecting human health or parasite transmission - merits a deeper dive.
Given the complexity of malaria, we needed a mathematical framework for
malaria epidemiology that could expose to scrutiny the relationship
between processes and patterns in human populations exposed to malaria
parasites as they age, with different intensities and with different
seasonal patterns. This software package, **`ramp.falciparum,`**
implements a new computational approach to malaria epidemiology (in that
narrow sense) using random variables. These methods are highly mimetic,
but computationally intensive and difficult to apply. To address these
limitations, we use hybrid variables to build a bridge from this
probabilistic approach to a set of simpler approximating models.

The epidemiology of *Plasmodium falciparum* malaria presents a unique
set of challenges due to the complex dynamics of infection, immunity,
disease and infectiousness as well as treatment and chemo-protection,
diagnostics and detection. Malaria can be measured in a dozen different
ways, but it has been difficult to present a simple synthesis of malaria
infection and disease in terms of the metrics that are commonly used in
research and clinical surveillance. An important metric is the
*Plasmodium falciparum* parasite rate, or *Pf*PR, defined as the average
prevalence of malaria taken from a cross-sectional survey. Another
metric, often measured as a covariate in research studies, is a parasite
count, the number of parasites in a blood slide field counted by a light
microscopist. In an old data set, collected during malariotherapy,
parasite counts fluctuated substantially over the time course of an
infection, but they were strongly statistically correlated with the *age
of infection* or AoI (Henry JM, *et al.*, 2022)[^1]. In malaria, the
*Pf*PR in several old studies had a characteristic shape when plotted
against age. Parasite densities have been used in research settings as
both a diagnostic criterion and as a correlate of disease. Malaria
epidemiology exhibits patterns that differ by diagnostic method, by
season, by sex, and by location.

With so many interacting factors, it has been a challenge to develop
model that could deal with everything all at once. One approach to
studying malaria infection has been to develop mechanistic models for
the dynamics of malaria infection within a single host. The most
prominent models of this type are OpenMalaria and eMod, but there have
been several others. These computational models made it possible to
develope comprehensive individual-based simulation models, or IBMs, for
malaria policy. While these approaches have been able to replicate the
patterns, the outputs of the models are usually just as complex as the
data collected in field studies. A synthesis of malaria epidemiology has
proven elusive.

[^1]: Henry JM, Carter A & Smith DL (2022) Infection age as a predictor
    of epidemiological metrics for malaria. Malar J 21; 117,
    <https://doi.org/10.1186/s12936-022-04134-5>

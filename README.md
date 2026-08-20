# `ramp.falciparum` <br><br> A Probabilistic Approach to Malaria Epidemiology

## Install

To install the latest version of `ramp.falciparum` from GitHub, run the following lines of code in an R session. We also recommend loading the supporting package `ramp.trace`

```
library(devtools)
devtools::install_github("dd-harp/ramp.falciparum")
devtools::install_github("dd-harp/ramp.trace")
```
Then load them: 
```
library(ramp.falciparum)
library(ramp.trace)
```

## Falciparum Malaria 

The health burden of falciparum malaria -- human malaria caused by *Plasmodium falciparum* -- is enormous, 
so there is a need to understand malaria epidemiology well enough to guide development of policies.
The epidemiology of falciparum malaria is also complex: the basic features of malaria include: 

+ Simple malaria infections have a complex **time course:** an acute growth phase (lasting up to a month) is followed by a chronic phase that lasts for months (on average). Over the time-course of infection, parasite densities fluctuate enormously. The toolbox for handling infections includes stage-of-infection, age-of-infection, within-host models, and probability theory. 

+ The rate of exposure to malaria often exceeds the slow rate that infections clear, so **superinfection** is common. *Queuing theory* has been the main toolbox: these are compartmental models that sub-divide infection into an infinite number of states -- the multiplicity of infection (MoI) -- and that explore the dynamics of the MoI. Here, we will also develop *hybrid variables* that model MoI distributions. 

+ Infections can be **treated and cured** with anti-malarial drugs, and treatment is often followed by a short period of **chemo-protection.** This raises questions about the adherence to the prescribed drug-regimens, and the evolution of drug resistance. 

+ **Infectiousness** -- the probability a mosquito becomes infected after blood feeding on a human -- is related to gametocyte densities and development of *transmission-blocking immunity.*  

+ **Disease** and the health burden of malaria is enormous; while there has been an inordinate focus on fever, the primary concerns are **

+ **Immunity** and its effects on infection, infectiousness, and disease; 

+ To understand malaria in populations, we must measure it. To be useful, it must be possible to relate the *true prevalence* of malaria to the prevalence that would be observed, the probability of **detection** by various **diagnostics.** 

+ If these models are to be used for policy, then it is essential to have outputs by **age.** 

All these are quantitative phenomena are interrelated, and there is probably no sensible way of describing them that does not involve mathematics.
Each one of these facets of malaria has been addressed in mathematical models, but it has proven difficult to formulate a synthesis.


This website holds code that implements mathematical models of malaria epidemiology.
The goal is to provide a single repository for mathematical models, and to develop some models of malaria that are useful for research and policy.
In some cases, the documentation includes some mathematical derivations.


## Robust Analytics 

Policy advice should be robust to uncertainty. To build rigor around the idea of robustness in policy, we recognized the need for a bespoke inferential framework designed for malaria analytics, which gave rise to **RAMP**
(Robust Analytics for Malaria Policy). An important goal of RAMP was to build an integrated computational environment to support a range of activities, including conventional statistical analysis and methods to characterize, quantify, and propagate uncertainty.
As RAMP took shape, RAMP software was being developed to implemented the principles in a stable computational form to facilitate the transformation of data into robust policy advice.
Today, core RAMP software includes several R packages (see [ramp.malaria](https://dd-harp.github.io/ramp.malaria/)) to deal with malaria epidemiology, transmission dynamics, and control in the broad sense. The [ramp.library](https://dd-harp.github.io/ramp.library/) includes reusable code for one of several models of human malaria infection and immuno-epidemiology. 

We recognized that the complex epidemiology of malaria - defined in the narrow sense to include only parasite infection in humans and processes affecting  human health or parasite transmission - merits a deeper dive.
Given the complexity of malaria, we needed a mathematical framework for malaria epidemiology that could expose to scrutiny the relationship between processes and patterns in human populations exposed to malaria parasites as they age, with different intensities and with different seasonal patterns.
This software package, **`ramp.falciparum,`** implements a new computational approach to malaria epidemiology (in that narrow sense) using random variables. These methods are highly mimetic, but computationally intensive and difficult to apply. To address these limitations, we use hybrid variables to build a bridge from this probabilistic approach to a set of simpler approximating models. 

The epidemiology of *Plasmodium falciparum* malaria presents a unique set of challenges due to the complex dynamics of infection, immunity, disease and infectiousness as well as treatment and chemo-protection, diagnostics and detection. Malaria can be measured in a dozen different ways, but it has been difficult to present a simple synthesis of malaria infection and disease in terms of the metrics that are commonly used in research and clinical surveillance. An important metric is the *Plasmodium falciparum* parasite rate, or *Pf*PR, defined as the average prevalence of malaria taken from a cross-sectional survey. Another metric, often measured as a covariate in research studies, is a parasite count, the number of parasites in a blood slide field counted by a light microscopist. In an old data set, collected during malariotherapy, parasite counts fluctuated substantially over the time course of an infection, but they were strongly statistically correlated with the *age of infection* or AoI (Henry JM, *et al.*, 2022)^[Henry JM, Carter A & Smith DL (2022) Infection age as a predictor of epidemiological metrics for malaria. Malar J 21; 117, https://doi.org/10.1186/s12936-022-04134-5]. In malaria, the *Pf*PR in several old studies had a characteristic shape when plotted against age. Parasite densities have been used in research settings as both a diagnostic criterion and as a correlate of disease. Malaria epidemiology exhibits patterns that differ by diagnostic method, by season, by sex, and by location.   

With so many interacting factors, it has been a challenge to develop model that could deal with everything all at once. One approach to studying malaria infection has been to develop mechanistic models for the dynamics of malaria infection within a single host. The most prominent models of this type are OpenMalaria and eMod, but there have been several others. These computational models made it possible to develope comprehensive individual-based simulation models, or IBMs, for malaria policy. While these approaches have been able to replicate the patterns, the outputs of the models are usually just as complex as the data collected in field studies. A synthesis of malaria epidemiology has proven elusive. 

We present the computational algorithms that support a probabilistic approach to malaria epidemiology.  We start with a semi-Markovian model of malaria exposure and infection, whose states are represented by random variables that describe the multiplicity of infection (MoI) in a host and the age of infection (AoI). Assuming that parasite densities can be predicted by the AoI in a statistical sense, we can compute probability distribution functions describing parasite densities, parasite counts, and detection in an individual chosen at random from the population. From this, we present a model for parasite detection and parasite counts. This same approach has been extended to predict disease, immunity, treatment with anti-malarial drugs, and a brief period of chemo-protection. 

The probabilistic approach is both highly realistic and descriptive, but our goal was a synthesis. This synthesis involves a few steps:

1. We develop formula and functions to compute the mean MoI, the mean AoI and all its moments, and the probability of detection. 

2. Hybrid models for the mean MoI for malaria superinfection were developed by Nåsell (1985)^[Nåsell I (1985). Hybrid Models of Tropical Infections, 1st edition. Springer-Verlag. https://doi.org/10.1007/978-3-662-01609-1], We extend this approach, developing systems of differential equations that track the mean and higher order moments of the distribution of the AoI. 

3. We a new random variable describing the age of the youngest infection (AoY). We show how the variable serves as a basis for computing parasite density distributions in complex infections. 

4. We derive a hybrid variable for the mean AoY. 

5. We demonstrate that a simple system of ordinary differential equations can be used in place of the random variables for most applications. 

To put it another way, we can reduce the behavior of these highly complex probabilistic systems to a simple system of equations that has a high degree of accuracy. The computational and conceptual simplicity of hybrid models have some simplicity over compartmental models and stochastic individual-based models, and with the supporting probabilistic framework, provide a sound basis for a synthesis of observational malaria epidemiology. 



## References



# `ramp.falciparum` <br><br> A Probabilistic Approach to Malaria Epidemiology 

Policy advice should be robust to uncertainty. 
Robust Analytics for Malaria Policy (RAMP) was developed as a bespoke inferential system for malaria that aims to characterize, quantify, and propagate uncertainty. 
RAMP software has implemented those principles in a stable computational form to facilitate the transformation of data into robust policy advice.
The core RAMP software includes four R packages, including several models of human malaria infection and immuno-epidemiology. We recognize that the complex epidemiology of malaria - defined narrowly to focus on parasite infection in humans and related within-human phenomena - merits a deeper dive.
Given the complexity of malaria, we needed a mathematical framework for malaria epidemiology that could expose to scrutiny the relationship between processes and patterns in populations.
This software package, `ramp.falciparum,` implements a new computational approach to malaria epidemiology, narrowly defined. Here, we describe malaria epidemiology using random variables, and we use hybrid variables to build a bridge from this probabilistic approach to a set of simpler approximating models. 

The epidemiology of *Plasmodium falciparum* malaria presents a unique set of challenges due to the complex dynamics of infection, immunity, disease and infectiousness as well as treatment and chemo-protection, diagnostics and detection. Malaria can be measured in a dozen different ways, but it has been difficult to present a simple synthesis of malaria infection and disease in terms of the metrics that are commonly used in research and clinical surveillance. An important metric is the *Plasmodium falciparum* parasite rate, or *Pf*PR, defined as the average prevalence of malaria taken from a cross-sectional survey. Another metric, often measured as a covariate in research studies, is a parasite count, the number of parasites in a blood slide field counted by a light microscopist. In an old data set, collected during malariotherapy, parasite counts fluctuated substantially over the time course of an infection, but they were strongly statistically correlated with the *age of infection* or AoI (Henry JM, *et al.*, 2022)^[Henry JM, Carter A & Smith DL (2022) Infection age as a predictor of epidemiological metrics for malaria. Malar J 21; 117, https://doi.org/10.1186/s12936-022-04134-5]. In malaria, the *Pf*PR in several old studies had a characterstic shape when plotted against age. Parasite densities have been used in research settings as both a diagnostic criterion and as a correlate of disease. Malaria epidemiology exhibits patterns that differ by diagnostic method, by season, by sex, and by location.   

With so many interacting factors, it has been a challenge to develop model that could deal with everything all at once. One approach to studying malaria infection has been to develop mechanistic models for the dynamics of malaria infection within a single host. The most prominent models of this type are OpenMalaria and eMod, but there have been several others. These computational models made it possible to develope comprehensive individual-based simulation models, or IBMs, for malaria policy. While these approaches have been able to replicate the patterns, the outputs of the models are usually just as complex as the data collected in field studies. A synthesis of malaria epidemiology has proven elusive. 

We present the computational algorithms that support a probabilistic approach to malaria epidemiology.  We start with a semi-Markovian model of malaria exposure and infection, whose states are represented by random variables that describe the multiplicity of infection (MoI) in a host and the age of infection (AoI). Assuming that parasite densities can be predicted by the AoI in a statistical sense, we can compute probability distribution functions describing parasite densities, parasite counts, and detection in an individual chosen at random from the population. From this, we present a model for parasite detection and parasite counts. This same approach has been extended to predict disease, immunity, treatment with anti-malarial drugs, and a brief period of chemo-protection. 

The probabilistic approach is both highly realistic and descriptive, but our goal was a synthesis. This synthesis involves a few steps:

1. We develop formula and functions to compute the mean MoI, the mean AoI and all its moments, and the probability of detection. 

2. Hybrid models for the mean MoI for malaria superinfection were developed by Nåsell (1985)^[Nåsell I (1985). Hybrid Models of Tropical Infections, 1st edition. Springer-Verlag. https://doi.org/10.1007/978-3-662-01609-1], We extend this approach, deeloping systems of differential equations that track the mean and higher order moments of the distribution of the AoI. 

3. We a new random variable describing the age of the youngest infection (AoY). We show how the variable serves as a basis for computing parasite density distributions in complex infections. 

4. We derive a hybrid variable for the mean AoY. 

5. We demonstrate that a simple system of ordinary differential equations can be used in place of the random variables for most applications. 

To put it another way, we can reduce the behavior of these highly complex probabilistic systems to a simple system of equations that has a high degree of accuracy. The computational and conceptual simplicity of hybrid models have some simplicity over compartmental models and stochastic individual-based models, and with the supporting probabilistic framework, provide a sound basis for a synthesis of observational malaria epidemiology. 

## Installation

To install and use the latest version from GitHub, run the following lines of code in an R session.

```
library(devtools)
devtools::install_github("dd-harp/ramp.falciparum")
install.packages("ramp.falciparum")
library(ramp.falciparum)
```


## References



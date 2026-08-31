# `ramp.falciparum` 

## Falciparum Malaria Mathematical Epidemiology

***

**INSTALL** 

To install the latest version of `ramp.falciparum` from GitHub, run the following lines of code in an R session. We also recommend loading the supporting package `ramp.func`

```
library(devtools)
devtools::install_github("dd-harp/ramp.falciparum")
devtools::install_github("dd-harp/ramp.func")
```
Then load them: 
```
library(ramp.falciparum)
library(ramp.func)
```

## Falciparum Malaria 

The burden of falciparum malaria -- human disease caused by infection with *Plasmodium falciparum* -- on human health and health systems is enormous. 
In countries where malaria is endemic, the prevalence of malaria varies enormously. 
In some places, most people test positive and more than half of all outpatients get diagnosed with malaria. 
Severe malaria is uncommon, but it can result in death. Mild disease -- subjective fever or malaise -- is quite common. 

There is a need to understand malaria epidemiology well enough to manage it, 
but the epidemiology of falciparum malaria is complex, 
and a synthesis has proven elusive. 
The management of malaria has two main goals: 
there is a need to implement policies that **reduce the burden** of disease; 
and most countries have the long term goal of **eliminating malaria.** 
To manage malaria, we must measure it. 
Malaria epidemiology is now more than a dozen dozen years old, if we start counting from Laveran. 
Even if we define malaria epidemiology narrowly to include only factors related to infection, 
excluding malaria transmission dynamics and control, we can identify these core issues: 

1. exposure;

2. the complex time course of infection and 

3. superinfection; 

4. heterogeneous exposure; 

5. disease, including fever, anemia, and severe disease;

6. the causes of death from severe disease; 

7. gametocytes and infectiousenss;

8. treatment with anti-malarial drugs and its effects on infection, infectiousness, and disease; 

9. immunity and its effects on infection, infectiousness, and disease; 

10. malaria in pregnancy;

11. vaccines

12. the problem of diagnostics and detection by light microscopy, rapid diagnostic tests (RDTs), and various alternatives; 

13. human age and sex and its effects on malaria; 

14. human genetic blood disorders associated with malaria and their effects on malaria; 

15. malaria spatial dynamics; 

16. other kinds of epidemiologically relevant heterogeneity, including nutritional status and its effects on malaria, coinfection with other pathogens, propensity to seek care, and adherance to antimalarial drug regimens;

17. genetic diversity in the parasite population as an underlying explanation for malaria epidemiology; and

18. the evolution of resistance to drugs, diagnostics, and vaccines.

Mathematical models of malaria can serve an indispensible role in all of this, in part, 
because they make it possible to convey quantitative information accurately. 
The models can also be useful for understanding malaria, and they have played an important role in malaria research and analytics. 
The first mathematical models of malaria were published in 1908 and 1911 by Ronald Ross, who wrote:

> *... all epidemiology, concerned as it is with the variation of disease from time to time or from place to place, must be considered mathematically, however many variables are implicated, if it is to be considered scientifically at all. To say that a disease depends upon certain factors is not to say much, until we can also form an estimate as to how largely each factor influences the whole result. And the mathematical method of treatment is really nothing but the application of careful reasoning to the problems at issue.*

In 1957, the preface to George Macdonald's book *The Epidemiology and Control of Malaria,*  he argues that the models had, so far, fallen short: 

> *The mathematical studies of Ross appeared an attractive approach to new explanation, but experiment showed they did not complete the picture or provide explanation...*

Macdonald noted that all of the work that had been done on the mathematical epidemiology of malaria --- by Ross, Lotka, and McKendrick --- had all been done on a single, simple model. 

Innovation in malaria mathematical epidemiology was slow, but it started to accelerate after the Global Malaria Eradication Programme (GMEP) ended in 1969.
In 1950, Macdonald had published a mathematical model of superinfection, but (as explained beautifully by Paul Fine), the mathematics were flawed.
Macdonald's superinfection math can be understood through a branch of mathematics called *queuing theory,* but even if he had got it right, 
it was one feature among many. 
The first model of malaria with immunity was published in 1974 as part of the Garki Project.
The first within-host model was published by Hellriegel in 1992. 
Severe disease first appeared in 1999.
Published models of malaria now consider most of the core issues listed above, and several individual-based models of malaria have been published,
but there is a need for synthesis. 

## Mathematical Epidemiology

This R package has the goal of making the mathematical epidemiology of malaria more accessible to everyone, 
and to explain some recent innovations that have made it possible to reduce the enormous complexity of malaria. 
In some ways, it is a very well developed sandbox. 

In developing this approach to malaria, we consider the ontogeny of malaria epidemiology. 
The mathematical formulations can be reduced to a simple question:

> **Given a history of exposure to malaria in a person of some age, what are the factors that determine the outcomes of the next infection?**

In taming the complexity, we use **random variables.** 
In forging a synthesis, we make extensive use of **hybrid variables** describing the distributions of 
the age of infection (AoI), the multiplicity of infection (MoI), and cumulative exposure.  

## Robust Analytics 

To be useful in policy the models must be 
fit for purpose, 
realistic enough to be compelling, 
accurate enough to useful, 
and simple enough to be insightful. 

If we want to use models effectively in malaria policy, we 
will need to go one step further and show that are models are realistic and accurate enough.
Some of the hard work can be done through model-model comparison.

These issues fall under the category of **robust analytics for malaria policy,**  that consider malaria epidemiology in a much broader sense, 
including malaria transmission dynamics and control. 
For malaria analytics, we have have developed some related software, including [**`ramp.xds`**](https://dd-harp.github.io/ramp.xds/) and a suite of R packages developed for robust, simulation-based analytics called [**SimBA**](https://faculty.washington.edu/smitdave/simba/) 


## References



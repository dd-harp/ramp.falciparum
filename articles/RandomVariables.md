# Random Variable Models

We present the computational algorithms that support a probabilistic
approach to malaria epidemiology. We start with a semi-Markovian model
of malaria exposure and infection, whose states are represented by
random variables that describe the multiplicity of infection (MoI) in a
host and the age of infection (AoI). Assuming that parasite densities
can be predicted by the AoI in a statistical sense, we can compute
probability distribution functions describing parasite densities,
parasite counts, and detection in an individual chosen at random from
the population. From this, we present a model for parasite detection and
parasite counts. This same approach has been extended to predict
disease, immunity, treatment with anti-malarial drugs, and a brief
period of chemo-protection.

The probabilistic approach is both highly realistic and descriptive, but
our goal was a synthesis. This synthesis involves a few steps:

1.  We develop formula and functions to compute the mean MoI, the mean
    AoI and all its moments, and the probability of detection.

2.  Hybrid models for the mean MoI for malaria superinfection were
    developed by Nåsell (1985)[^1], We extend this approach, developing
    systems of differential equations that track the mean and higher
    order moments of the distribution of the AoI.

3.  We a new random variable describing the age of the youngest
    infection (AoY). We show how the variable serves as a basis for
    computing parasite density distributions in complex infections.

4.  We derive a hybrid variable for the mean AoY.

5.  We demonstrate that a simple system of ordinary differential
    equations can be used in place of the random variables for most
    applications.

To put it another way, we can reduce the behavior of these highly
complex probabilistic systems to a simple system of equations that has a
high degree of accuracy. The computational and conceptual simplicity of
hybrid models have some simplicity over compartmental models and
stochastic individual-based models, and with the supporting
probabilistic framework, provide a sound basis for a synthesis of
observational malaria epidemiology.

[^1]: Nåsell I (1985). Hybrid Models of Tropical Infections, 1st
    edition. Springer-Verlag.
    <https://doi.org/10.1007/978-3-662-01609-1>

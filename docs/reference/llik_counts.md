# Compute the log likelihood for a set of parasites counts data

Compute the log likelihood for a set of parasites counts data

## Usage

``` r
llik_counts(prms, alpha, counts, par_Fmu, par_O, F_b, F_Sa, q = 6, bvm = 13)
```

## Arguments

- prms:

  A set of parameters to be fitted

- alpha:

  the AoI for the observations

- counts:

  the counts

- par_Fmu:

  an object to dispatch
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- par_O:

  an object to dispath
  [d_Omega](https://dd-harp.github.io/ramp.falciparum/reference/d_Omega.md)

- F_b:

  a function to constrain the upper bound for fitting for
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- F_Sa:

  a function to constrain the slope for fitting for
  [Fmu](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.md)

- q:

  the sampling volume (in \\\log\_{10} microliters\\)

- bvm:

  the blood volume (in \\\log\_{10} microliters\\)

## Value

the negative log likelihood for a set of observations

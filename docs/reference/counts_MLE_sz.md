# Compute the MLE for parasite counts data

Compute the MLE for parasite counts data

## Usage

``` r
counts_MLE_sz(
  alpha,
  counts,
  inits,
  sz,
  F_b,
  F_Sa,
  q = 6,
  bvm = 13,
  SANN = FALSE
)
```

## Arguments

- alpha:

  the ages of infection

- counts:

  the counts

- inits:

  the initial guesses for the parameters

- sz:

  the size parameter

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

- SANN:

  if true, then use simulated annealing

## Value

the fitted parameter values and the log likelihood

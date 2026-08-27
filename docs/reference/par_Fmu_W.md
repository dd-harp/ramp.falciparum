# Set up parameters for [Fmu.W](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.W.md)

Set up parameters for
[Fmu.W](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.W.md)

## Usage

``` r
par_Fmu_W(
  peak = 20,
  liver = 7,
  tildeb = 10.3,
  tildel = 2,
  Sa = 0.0033,
  Sw = 0.001
)
```

## Arguments

- peak:

  The age of infection (in days) when parasite densities peak

- liver:

  The age of infection (in days) when parasites emerge from the liver

- tildeb:

  The maximum expected log10 parasite densities

- tildel:

  The minimum expected log10 parasite densities

- Sa:

  The decline in mu with respect to alpha

- Sw:

  The decline in mu with respect to immunity

## Value

mean log10 parasite densities

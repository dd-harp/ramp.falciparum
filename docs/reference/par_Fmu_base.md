# Set up parameters for [Fmu.base](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.base.md)

Set up parameters for
[Fmu.base](https://dd-harp.github.io/ramp.falciparum/reference/Fmu.base.md)

## Usage

``` r
par_Fmu_base(peak = 20, liver = 7, tildeb = 10, tildel = 2, Sa = 0.003)
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

## Value

a [list](https://rdrr.io/r/base/list.html)

# The random generation function for the age of the youngest infection (AoY)

The random generation function for the age of the youngest infection
(AoY)

## Usage

``` r
rAoY(N, a, FoI_a, bday = 0, hhat = 1, r = 1/200, alphamin = 0)
```

## Arguments

- N:

  the number of observations

- a:

  cohort age

- FoI_a:

  a cohort trace function

- bday:

  the cohort birthday

- hhat:

  a local scaling parameter for the FoI

- r:

  the clearance rate for a simple infection

- alphamin:

  the minimum value of the AoI to return

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(N)

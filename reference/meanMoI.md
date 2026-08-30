# The mean MoI in a host cohort of age \\a\\

The mean multiplicity of infection (MoI) is \$\$m\_\bday(a\|h) =
\int_0^a z\_\bday(\alpha, a\|h) d \alpha\$\$

## Usage

``` r
meanMoI(a, FoI_a, bday = 0, hhat = 1, r = 1/200)
```

## Arguments

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

## Value

a [numeric](https://rdrr.io/r/base/numeric.html) vector of length(a)

# Potter Cove Body Mass and Biomass Dataset

This dataset contains curated body mass, biomass, and resource density
information for consumers and resources in the Potter Cove food web,
from Rodriguez & Saravia (2024). The dataset was processed from raw
measurements, with units converted to kilograms (kg), and some special
handling for detritus and sponge taxa where direct body mass or density
values were not available or meaningful.

## Usage

``` r
PotterCove_bm
```

## Format

A data frame with the following columns:

- consumer:

  Name of the consumer species.

- resource:

  Name of the resource species.

- c_bodymass:

  Body mass of the consumer (in kg).

- r_bodymass:

  Body mass of the resource (in kg). May be set to `1` for detritus or
  sponges.

- r_biomass:

  Biomass of the resource (in kg).

- r_density:

  Resource density (units depend on context). May be `1` for some
  resources or `-999` to flag missing values.

## Source

<https://doi.org/10.5281/zenodo.10790590>

## Details

When resource density is unknown, the value `-999` is used as a flag so
that downstream functions can infer density or body mass using
allometric relationships following Pawar et al. (2012). Body mass and
biomass were originally in grams and converted to kilograms. Sponge taxa
(e.g., *Chalinidae*, *Dendrilla antarctica*) resource body mass and
density were assigned to 1.

## References

Rodriguez, I. D., & Saravia, L. A. (2024). Potter Cove’s Heavyweights:
Estimation of Species’ Interaction Strength of an Antarctic Food Web.
Ecology and Evolution, 14(11), e70389.
https://doi.org/10.1002/ece3.70389

Pawar, S., Dell, A.I., & Savage, V.M. (2012). Dimensionality of consumer
search space drives trophic interaction strengths. *Nature*, 486(7404),
485–489.

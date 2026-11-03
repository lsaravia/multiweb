# Generate Node aligned Interlayer Edges for a Multilayer Network from the function `convert_to_intra_format`

Generate Node aligned Interlayer Edges for a Multilayer Network from the
function `convert_to_intra_format`

## Usage

``` r
add_interlayer_links(intra_data, interlayer_weight = 1)
```

## Arguments

- intra_data:

  A list with `$vertices` (node info) and `$intra` (intralayer edges).

- interlayer_weight:

  Numeric, the weight for interlayer links (default: 1).

## Value

A dataframe with interlayer edges.

# Compute Centrality for a Weighted/Unweighted Network

This function calculates the centrality of species in a network based on
an adjacency matrix or igraph object using a specified centrality
measure.

## Usage

``` r
calc_centrality(A, centrality_func = eigen_centrality)
```

## Arguments

- A:

  An igraph object or an adjacency matrix. The values represent
  interaction strengths.

- centrality_func:

  A function to compute centrality. Default is `eigen_centrality`. Other
  options: `page_rank`, `betweenness`, `degree`, etc.

## Value

A list containing:

- CentralityScores:

  A data frame with species names and their centrality scores.

- Plot_Centrality:

  A ggplot object showing the top 10 species by centrality.

## Examples

``` r
# Example usage with an adjacency matrix A
results <- calc_centrality(A, centrality_func = page_rank)
#> Error: object 'A' not found
print(results$CentralityScores)
#> Error: object 'results' not found
print(results$Plot_Centrality)
#> Error: object 'results' not found
```

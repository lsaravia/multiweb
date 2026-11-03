# From Generalized Lotka Volterra adjacency matrix to igraph topological object

It counts predator-prey/Antagonistic interactions like 1 edge,
competition and mutualisms are counted as is: two edges or one edge

## Usage

``` r
fromGLVadjToIgraph(glvAdj, spc)
```

## Arguments

- glvAdj:

  numeric matrix Generalized Lotka-Volterra adjacency matrix

- spc:

  numeric vector of present Species

## Value

an igraph object

## Examples

``` r
# Build a matrix

m <- matrix(0,nrow=5,ncol=5)
m[1,2] <- m[1,3] <- m[3,4]<- .2
m[2,1] <- m[3,1] <- m[4,3] <- -2
m[5,4] <- m[4,5] <- 0.1            # Mutualistic
m[1,1] <- -0.01                    # Cannibalistic

g <- fromGLVadjToIgraph(m,c(1,1,1,1,0))

g <- fromGLVadjToIgraph(m,c(0,1,1,1,1))
```

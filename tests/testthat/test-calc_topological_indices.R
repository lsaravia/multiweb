
test_that("calc_topological_indices", {
    library(igraph)
    g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, 2-+7, simplify = FALSE)
    expect_snapshot(
     calc_topological_indices(g)
    )
})

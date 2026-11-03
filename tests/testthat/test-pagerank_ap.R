test_that("calc_pagerank_ap ranks resources above consumers in a chain", {
  g <- igraph::graph_from_literal(Resource -+ Herbivore, Herbivore -+ Predator)
  r <- calc_pagerank_ap(g)
  expect_equal(r$Species, c("Resource", "Herbivore", "Predator"))
  expect_equal(r$PageRank, c(1/2, 1/3, 1/6), tolerance = 1e-8)
  expect_equal(sum(r$PageRank), 1)
})

test_that("calc_pagerank_ap matches igraph PageRank on the reversed graph with root", {
  g <- netData[[1]]
  r <- calc_pagerank_ap(g)
  S <- igraph::vcount(g)
  ga <- igraph::add_vertices(g, 1, name = "ROOT")
  b <- which(igraph::degree(g, mode = "in") == 0)
  ga <- igraph::add_edges(ga, as.vector(rbind(S + 1, b)))
  ga <- igraph::add_edges(ga, as.vector(rbind(1:S, S + 1)))
  ga <- igraph::simplify(ga)
  p <- suppressWarnings(igraph::page_rank(igraph::reverse_edges(ga), damping = 0.999999)$vector[1:S])
  p <- p / sum(p)
  expect_equal(unname(p[r$Species]), r$PageRank, tolerance = 1e-5)
})

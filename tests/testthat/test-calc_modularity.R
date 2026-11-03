library(testthat)
library(igraph)
library(dplyr)

test_that("calc_modularity returns correct output format for single network", {
  # Create an undirected and unweighted test graph
  g <- netData[[1]]

  # Run function
  result <- calc_modularity(g, cluster_function = cluster_spinglass)

  # Check output format
  expect_s3_class(result, "data.frame")
  expect_true("Modularity" %in% colnames(result))
  expect_type(result$Modularity, "double")
})

test_that("calc_modularity handles directed and weighted networks", {
  # Create a directed and weighted test graph
  g <- netData[[1]]
  E(g)$weight <- runif(ecount(g))  # Assign random edge weights

  # Run function with weights
  result <- calc_modularity(g, cluster_function = cluster_infomap)

  # Check output format
  expect_s3_class(result, "data.frame")
  expect_true("Modularity" %in% colnames(result))
  expect_type(result$Modularity, "double")
})

test_that("calc_modularity propagates errors from clustering functions that do not support directed graphs", {
  g <- netData[[1]]
  # cluster_louvain only works on undirected graphs
  expect_error(calc_modularity(g, cluster_function = cluster_louvain),
               "undirected")
})

test_that("calc_modularity computes modularity in the correct range", {
  g <- netData[[1]]

  # Run function
  result <- calc_modularity(g, cluster_function = cluster_spinglass)

  # Modularity of the spinglass partition recomputed on the directed network (~0.19-0.21)
  expect_true(all(result$Modularity >= 0.18 & result$Modularity <= 0.22))
})

test_that("calc_modularity works with both cluster_spinglass and cluster_infomap", {
  g <- netData[[1]]

  result1 <- calc_modularity(g, cluster_function = cluster_spinglass)
  expect_s3_class(result1, "data.frame")

  result2 <- calc_modularity(g, cluster_function = cluster_infomap)
  expect_s3_class(result2, "data.frame")
})

test_that("calc_modularity handles a list of networks", {
  skip_on_cran()  # Avoid running on CRAN due to randomness

  # Generate a list of randomized networks from netData
  nullg <- generateERbasal(netData[[1]], 10)

  # Run function on list
  result <- calc_modularity(nullg, cluster_function = cluster_spinglass)

  # Check output format
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == length(nullg))  # Should return one row per network
})

test_that("calc_modularity handles parallel processing", {
  g <- generateERbasal(netData[[1]], 10)

  # Run sequential
  result_seq <- calc_modularity(g, ncores = 0)
  expect_s3_class(result_seq, "data.frame")

  # Run parallel
  result_par <- calc_modularity(g, ncores = 2)
  expect_s3_class(result_par, "data.frame")
})

# Two 8-node cliques joined by a few cross edges. With heavy intra-clique weights
# and weak inter-clique weights both algorithms recover the planted partition.
make_planted_graph <- function() {
  n <- 8
  within <- function(off) t(combn(seq_len(n) + off, 2))
  cross <- cbind(c(1, 3, 5, 7), c(9, 12, 14, 16))
  el <- rbind(within(0), within(n), cross)
  g <- graph_from_edgelist(el, directed = FALSE)
  V(g)$name <- paste0("sp", seq_len(vcount(g)))
  planted <- rep(1:2, each = n)
  is_cross <- planted[el[, 1]] != planted[el[, 2]]
  list(g = g, planted = planted, w = ifelse(is_cross, 0.1, 10))
}

test_that("calc_modularity uses the 'weight' edge attribute with cluster_infomap", {
  p <- make_planted_graph()
  gw <- p$g
  E(gw)$weight <- p$w

  expected_w  <- modularity(p$g, p$planted, weights = p$w)
  expected_uw <- modularity(p$g, p$planted)
  expect_gt(abs(expected_w - expected_uw), 0.01)  # the test can tell them apart

  expect_equal(calc_modularity(gw, cluster_function = cluster_infomap)$Modularity,
               expected_w, tolerance = 1e-6)
  expect_equal(calc_modularity(p$g, cluster_function = cluster_infomap)$Modularity,
               expected_uw, tolerance = 1e-6)
})

test_that("calc_modularity uses the 'weight' edge attribute with cluster_spinglass", {
  p <- make_planted_graph()
  gw <- p$g
  E(gw)$weight <- p$w

  expect_equal(calc_modularity(gw, cluster_function = cluster_spinglass)$Modularity,
               modularity(p$g, p$planted, weights = p$w), tolerance = 1e-6)
})

test_that("calc_modularity returns one weighted value per network in a list", {
  p <- make_planted_graph()
  gw <- p$g
  E(gw)$weight <- p$w

  res <- calc_modularity(list(gw, p$g), cluster_function = cluster_infomap)
  expect_equal(nrow(res), 2)
  expect_equal(res$Modularity,
               c(modularity(p$g, p$planted, weights = p$w),
                 modularity(p$g, p$planted)),
               tolerance = 1e-6)
})

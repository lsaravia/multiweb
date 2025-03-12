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
  result <- calc_modularity(g, weights = E(g)$weight, cluster_function = cluster_infomap)

  # Check output format
  expect_s3_class(result, "data.frame")
  expect_true("Modularity" %in% colnames(result))
  expect_type(result$Modularity, "double")
})

test_that("calc_modularity only accepts cluster_spinglass and cluster_infomap", {
  g <- netData[[1]]
  # Expect error for an invalid clustering function
  expect_error(calc_modularity(g, cluster_function = cluster_louvain),
               "must be either 'igraph::cluster_spinglass' or 'igraph::cluster_infomap'")
})

test_that("calc_modularity computes modularity in the correct range", {
  g <- netData[[1]]

  # Run function
  result <- calc_modularity(g, cluster_function = cluster_spinglass)

  # Modularity values should be in range (-0.5 to 1)
  expect_true(all(result$Modularity >= 0.18 & result$Modularity <= 0.19))
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

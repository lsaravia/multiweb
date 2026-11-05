test_that("curve_ball preserves in and out degree and returns igraph objects", {
  g <- netData[[1]]
  set.seed(1)
  nets <- curve_ball(g, nsim = 5)
  expect_length(nets, 5)
  for (x in nets) {
    expect_s3_class(x, "igraph")
    expect_equal(unname(degree(x, mode = "out")), unname(degree(g, mode = "out")))
    expect_equal(unname(degree(x, mode = "in")), unname(degree(g, mode = "in")))
    expect_equal(components(x)$no, 1)  # the original has 1 component
  }
})

test_that("curve_ball really randomizes the network", {
  g <- netData[[1]]
  A0 <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  set.seed(1)
  kept <- sapply(curve_ball(g, nsim = 20), function(x) {
    A <- as.matrix(as_adjacency_matrix(x, sparse = FALSE))
    sum(A & A0) / sum(A0)
  })
  # With 5 trades (bug) ~93% of the links were kept; fully randomized is ~40%
  expect_lt(mean(kept), 0.6)
})

test_that("curve_ball does not avoid the diagonal (cannibalism is allowed)", {
  g <- netData[[1]]
  expect_gt(sum(which_loop(g)), 0)
  set.seed(1)
  loops <- sapply(curve_ball(g, nsim = 30), function(x) sum(which_loop(x)))
  expect_gt(max(loops), 0)
})

test_that("curve_ball and shuffle_network_deg give similar degree-preserving randomizations", {
  g <- netData[[1]]
  A0 <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  set.seed(1)
  kept_cb <- sapply(curve_ball(g, nsim = 20), function(x) {
    A <- as.matrix(as_adjacency_matrix(x, sparse = FALSE)); sum(A & A0) / sum(A0)
  })
  kept_sh <- sapply(1:20, function(i) {
    A <- shuffle_network_deg(g, delta = 20 * ecount(g), weighted = FALSE)
    sum((A > 0) & A0) / sum(A0)
  })
  expect_lt(abs(mean(kept_cb) - mean(kept_sh)), 0.05)
})

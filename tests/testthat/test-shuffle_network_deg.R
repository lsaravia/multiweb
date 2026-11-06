# Random directed network with `nloops` self-loops
make_loop_net <- function(n = 40, conn = 0.15, nloops = 8) {
  A <- matrix(0, n, n)
  off <- which(row(A) != col(A))
  A[sample(off, round(conn * n * n))] <- 1
  diag(A)[sample(n, nloops)] <- 1
  A
}

# Copy of shuffle_network_deg() before the fix_loops argument (unweighted case)
shuffle_old <- function(A, delta) {
  for (iter in 1:delta) {
    edges <- which(A > 0, arr.ind = TRUE)
    if (nrow(edges) < 2) break
    s <- sample(1:nrow(edges), 2)
    i <- edges[s[1], 1]; j <- edges[s[1], 2]
    i_prime <- edges[s[2], 1]; j_prime <- edges[s[2], 2]
    if (i != i_prime && j != j_prime && A[i, j_prime] == 0 && A[i_prime, j] == 0) {
      A[i, j_prime] <- A[i, j]
      A[i_prime, j] <- A[i_prime, j_prime]
      A[i, j] <- 0
      A[i_prime, j_prime] <- 0
    }
  }
  A
}

test_that("fix_loops = TRUE keeps the diagonal fixed", {
  set.seed(1)
  A <- make_loop_net()
  B <- shuffle_network_deg(A, delta = 500, weighted = FALSE, fix_loops = TRUE)
  expect_identical(diag(B), diag(A))
  expect_false(isTRUE(all.equal(B, A)))  # but the rest was randomized
})

test_that("fix_loops = FALSE can change the loops", {
  set.seed(1)
  A <- make_loop_net()
  loops <- sapply(1:20, function(i) sum(diag(shuffle_network_deg(A, delta = 500, weighted = FALSE))))
  expect_true(any(loops != sum(diag(A))))
})

test_that("degrees and number of links are preserved, without duplicated links", {
  set.seed(2)
  A <- make_loop_net()
  for (fl in c(TRUE, FALSE)) {
    B <- shuffle_network_deg(A, delta = 500, weighted = FALSE, fix_loops = fl)
    expect_equal(rowSums(B), rowSums(A))
    expect_equal(colSums(B), colSums(A))
    expect_equal(sum(B), sum(A))
    expect_true(all(B %in% c(0, 1)))  # no cell with more than one link
  }
})

test_that("fix_loops = FALSE reproduces the previous version with the same seed", {
  set.seed(3)
  A <- make_loop_net()
  set.seed(10); old <- shuffle_old(A, 300)
  set.seed(10); new <- shuffle_network_deg(A, delta = 300, weighted = FALSE)
  expect_identical(new, old)
})

test_that("fix_loops works with igraph input and weights", {
  set.seed(4)
  A <- make_loop_net(n = 20, conn = 0.2, nloops = 4)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  igraph::E(g)$weight <- runif(igraph::ecount(g), 1, 10)
  W <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  B <- shuffle_network_deg(g, delta = 200, weighted = TRUE, fix_loops = TRUE)
  expect_equal(diag(B), diag(W))
  expect_equal(rowSums(B > 0), rowSums(W > 0))
  expect_equal(colSums(B > 0), colSums(W > 0))
  expect_equal(sort(B[B > 0]), sort(W[W > 0]))  # same weights
})

test_that("shuffle_network_deg warns for a multigraph when unweighted", {
  A <- make_loop_net(n = 10, conn = 0.3, nloops = 0)
  A[which(A > 0)[1]] <- 2
  expect_warning(shuffle_network_deg(A, delta = 5, weighted = FALSE), "multigraph")
  expect_no_warning(shuffle_network_deg(A, delta = 5, weighted = TRUE))  # values are weights
})

test_that("curve_ball fix_loops keeps loops fixed and preserves degrees", {
  g <- netData[[1]]
  A0 <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  set.seed(5)
  nets <- curve_ball(g, nsim = 10, fix_loops = TRUE)
  for (x in nets) {
    A <- as.matrix(igraph::as_adjacency_matrix(x, sparse = FALSE))
    expect_identical(unname(diag(A)), unname(diag(A0)))
    expect_equal(unname(rowSums(A)), unname(rowSums(A0)))
    expect_equal(unname(colSums(A)), unname(colSums(A0)))
  }
})

test_that("curve_ball fix_loops with several components and many loops", {
  set.seed(6)
  A <- make_loop_net(n = 30, conn = 0.08, nloops = 6)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  nets <- curve_ball(g, nsim = 10, fix_loops = TRUE)
  for (x in nets) {
    B <- as.matrix(igraph::as_adjacency_matrix(x, sparse = FALSE))
    expect_identical(unname(diag(B)), unname(diag(A)))
    expect_equal(unname(rowSums(B)), unname(rowSums(A)))
    expect_equal(unname(colSums(B)), unname(colSums(A)))
  }
})

test_that("curve_ball warns for a multigraph", {
  g <- igraph::make_graph(c(1, 2, 1, 2, 2, 3, 3, 1), directed = TRUE)
  expect_warning(curve_ball(g, nsim = 1), "multigraph")
})

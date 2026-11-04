skip_if(Sys.which("infomap") == "", "infomap binary not installed")

# Two 8-node cliques joined by a few cross edges (see test-calc_modularity.R)
make_planted_graph <- function() {
  n <- 8
  within <- function(off) t(combn(seq_len(n) + off, 2))
  cross <- cbind(c(1, 3, 5, 7), c(9, 12, 14, 16))
  g <- graph_from_edgelist(rbind(within(0), within(n), cross), directed = FALSE)
  V(g)$name <- paste0("sp", seq_len(vcount(g)))
  g
}

new_outdir <- function() {
  d <- tempfile("infomap_")
  dir.create(d)
  d
}

test_that("run_infomap returns an igraph communities object", {
  g <- netData[[1]]
  res <- run_infomap(g, output_dir = new_outdir())

  expect_s3_class(res, "communities")
  expect_equal(length(membership(res)), vcount(g))
  expect_equal(names(membership(res)), V(g)$name)
  expect_true(is.numeric(res$codelength) && res$codelength > 0)
  expect_gte(modularity(res), -0.5)
  expect_lte(modularity(res), 1)
})

test_that("run_infomap return_df gives module, node and flow", {
  g <- netData[[1]]
  res <- run_infomap(g, output_dir = new_outdir(), return_df = TRUE)

  expect_named(res, c("communities", "codelength"))
  expect_named(res$communities, c("module", "node", "flow"))
  expect_equal(nrow(res$communities), vcount(g))
  expect_setequal(res$communities$node, V(g)$name)
  expect_equal(sum(res$communities$flow), 1, tolerance = 1e-3)
})

test_that("run_infomap recovers planted modules", {
  g <- make_planted_graph()
  res <- run_infomap(g, output_dir = new_outdir(), directed = FALSE)

  expect_equal(compare(membership(res), rep(1:2, each = 8), method = "nmi"), 1)
  expect_gt(modularity(res), 0.4)
})

test_that("run_infomap uses the 'weight' edge attribute", {
  g <- netData[[1]]
  gw <- g
  set.seed(2)
  E(gw)$weight <- runif(ecount(g), 0.01, 100)

  res_uw <- run_infomap(g,  output_dir = new_outdir())
  res_w  <- run_infomap(gw, output_dir = new_outdir())

  expect_false(isTRUE(all.equal(res_uw$codelength, res_w$codelength)))
})

test_that("run_infomap weights move a bridging node into the strongly linked module", {
  g <- make_planted_graph()
  # Node 1 (clique A) gets one very heavy edge to node 9 (clique B) and light ones inside A
  E(g)$weight <- 1
  for (j in 2:8) E(g, P = c(1, j))$weight <- 0.01
  E(g, P = c(1, 9))$weight <- 100

  res <- run_infomap(g, output_dir = new_outdir(), directed = FALSE)
  mem <- membership(res)
  expect_equal(unname(mem["sp1"]), unname(mem["sp9"]))
})

test_that("run_infomap is reproducible for a given seed", {
  g <- netData[[1]]
  r1 <- run_infomap(g, output_dir = new_outdir(), seed = 42)
  r2 <- run_infomap(g, output_dir = new_outdir(), seed = 42)
  expect_equal(membership(r1), membership(r2))
  expect_equal(r1$codelength, r2$codelength)
})

test_that("run_infomap stores the weighted modularity of its partition", {
  g <- netData[[1]]
  set.seed(2)
  E(g)$weight <- runif(ecount(g), 0.1, 10)

  res <- run_infomap(g, output_dir = new_outdir())
  mem <- membership(res)
  expect_equal(modularity(res), modularity(g, mem, weights = E(g)$weight, directed = TRUE))
  # it is not the unweighted value
  expect_false(isTRUE(all.equal(modularity(res), modularity(g, mem, directed = TRUE))))

  res_und <- run_infomap(g, output_dir = new_outdir(), directed = FALSE)
  expect_equal(modularity(res_und),
               modularity(g, membership(res_und), weights = E(g)$weight, directed = FALSE))
})

test_that("run_infomap stored modularity matches the unweighted partition modularity without weights", {
  g <- netData[[1]]
  res <- run_infomap(g, output_dir = new_outdir())
  expect_equal(modularity(res), modularity(g, membership(res), directed = TRUE))
})

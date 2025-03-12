

test_that("calc_topological_roles", {
  g <- netData[[1]]
  set.seed(1432)
  comm <- cluster_spinglass(g)
  expect_snapshot(
    calc_topological_roles(g,nsim=1,community=comm)
  )
})

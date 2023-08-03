

test_that("calc_topological_roles", {
  g <- netData[[1]]
  set.seed(1432)
  expect_snapshot(
    calc_topological_roles(g,nsim=10,ncores=4)
  )
})

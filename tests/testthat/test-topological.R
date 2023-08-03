test_that("calc and classify topological roles", {

  g <- netData[[1]]
  set.seed(1432)
  tp <- calc_topological_roles(g,nsim=10,ncores=4)
  expect_snapshot(
    classify_topological_roles(tp,g,plt=FALSE))
})

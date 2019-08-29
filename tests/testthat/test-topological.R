context("topological roles")

test_that("calc and classify topological roles", {

  g <- netData[[1]]

  tp <- calc_topological_roles(g,nsim=10,ncores=4)
  tmp <- tempfile()
  expect_known_output(classify_topological_roles(tp,g,plt=FALSE), tmp)

})

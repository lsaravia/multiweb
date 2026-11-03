test_that("calc and classify topological roles", {

  g <- netData[[1]]
  set.seed(1432)
  tp <- calc_topological_roles(g,nsim=10,ncores=4)
  res <- classify_topological_roles(tp,g)
  expect_named(res, c("hub_conn", "gp"))
  expect_s3_class(res$gp, "ggplot")
  expect_true(all(c("node", "node_name", "type") %in% names(res$hub_conn)))
  expect_true(all(res$hub_conn$type %in% c("hubcon", "modhub", "modcon", "modspe")))
})

# Test ReadNetwork
context("Read netwokr files")
library(EcoNetwork)

test_that("Read a single file", {
  fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
  expect_known_value(g <- readNetwork(fileName))

})

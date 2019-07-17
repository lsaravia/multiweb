# Test ReadNetwork
context("Read network files")
library(multiweb)

test_that("Read a single csv file with adyacency matrix format and header and first column names", {
  fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
  g <- readNetwork(fileName)
  expect_is(g,"igraph")
  expect_equal(igraph::vcount(g),180)
  expect_equal(igraph::ecount(g),1546)
})

test_that("Read a single csv file with interaction list format", {
  fileName <- system.file("extdata", "WeddellSea_FW.csv", package = "multiweb")
  g <- readNetwork(fileName)
  expect_is(g,"igraph")
  expect_equal(igraph::vcount(g),442)
  expect_equal(igraph::ecount(g),1915)
})

test_that("Read multiple files with compatible formats", {
  filepath <- c(system.file("extdata",  package = "multiweb"))

  dn <- list.files(filepath,pattern = "^.*\\.csv$")
  g <- readNetwork(dn,filepath)
  expect_is(g,"list")
  expect_equal(length(g),2)
  expect_match(names(g)[2],"WeddellSea_FW")
  expect_is(g[[2]],"igraph")
})


test_that("Read a single txt file adyacency matrix no header no species names", {
  fileName <- system.file("extdata", "carpinteria_FW.txt", package = "multiweb")
  # The file has no header but also no species names so a warning is expected
  expect_warning(readNetwork(fileName,fhead = FALSE))
  # No header and no species names
  g <-   readNetwork(fileName,fhead = FALSE,skipColumn = 0)
  expect_is(g,"igraph")
  expect_equal(igraph::vcount(g),273)
  expect_equal(igraph::ecount(g),3878)
})

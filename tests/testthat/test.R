test_that("test data", {
  data("CedarCreek.2004.2006.dat")
  expect_that(CedarCreek.2004.2006.dat, is_a("data.frame"))
}
)

test_that("test data", {
  data("CedarCreek.2004.2006.res")
  expect_that(CedarCreek.2004.2006.res, is_a("list"))
}
)

test_that("test data", {
  data("CedarCreek.2004.res")
  expect_that(CedarCreek.2004.2006.res, is_a("list"))
}
)

test_that("test data", {
  data("CedarCreek.2004.2006.dat")
  expect_that(CedarCreek.2004.2006.dat, is_a("data.frame"))
}
)

test_that("test res year 2004", {
  data("CedarCreek.2004.res")
  expect_that(CedarCreek.2004.res, is_a("list"))

  expect_equal(length(CedarCreek.2004.res$fobs),
               length(CedarCreek.2004.res$xpr))
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$mOccur)[1])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$mCal)[2])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$mPrd)[2])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$mMotifs)[2])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$tCal)[2])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$tPrd)[2])
  expect_equal(length(CedarCreek.2004.res$fobs),
               dim(CedarCreek.2004.res$tNbcl)[2])

  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$tree.I$aff)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$tree.I$aff)[2])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$mCal)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$mPrd)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$mMotifs)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$tCal)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$tPrd)[1])
  expect_equal(dim(CedarCreek.2004.res$mOccur)[2],
               dim(CedarCreek.2004.res$tNbcl)[1])
}
)

test_that("test res years 2004 to 2006", {

  data("CedarCreek.2004.2006.res")

  expect_that(CedarCreek.2004.2006.res, is_a("list"))

  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               length(CedarCreek.2004.2006.res$xpr))
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$mOccur)[1])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$mCal)[2])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$mPrd)[2])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$mMotifs)[2])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$tCal)[2])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$tPrd)[2])
  expect_equal(length(CedarCreek.2004.2006.res$fobs),
               dim(CedarCreek.2004.2006.res$tNbcl)[2])

  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$tree.I$aff)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$tree.I$aff)[2])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$mCal)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$mPrd)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$mMotifs)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$tCal)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$tPrd)[1])
  expect_equal(dim(CedarCreek.2004.2006.res$mOccur)[2],
               dim(CedarCreek.2004.2006.res$tNbcl)[1])
}
)

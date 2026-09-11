# `testthat` is the R package used to run the package's automated tests. It
# provides helpers such as `test_check()`, `test_that()`, and `expect_equal()`
# so the tests can state expected behavior in readable R code.
library(testthat)

# Load this package before running its tests, so the tests call the installed
# package functions in the same way a user would.
library(ewgroup)

# Ask `testthat` to find and run all files in tests/testthat/.
test_check("ewgroup")

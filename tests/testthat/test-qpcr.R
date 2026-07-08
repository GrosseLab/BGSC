context("comparativeMethod_qPCR.analysis")

test_that("relative expression is higher for a lower gene Ct (more starting template)", {
  Ref <- c(20, 20, 20)
  Gene <- c(25, 23, 27) # sample 1 is the calibrator (t0); sample 2 has a lower Ct -> more expression

  res <- comparativeMethod_qPCR.analysis(Gene, Ref)

  # ddCt[2] = (23-20) - (25-20) = -2  ->  fold change = 2^(-(-2)) = 4  ->  100 * 4 = 400
  expect_equal(res$relative_values[2], 400, tolerance = 1e-8)
  expect_equal(res$ddCTRelExp[2], 400, tolerance = 1e-8)

  # sample 3 has a higher Ct than the calibrator -> lower relative expression
  expect_lt(res$relative_values[3], 100)
  expect_lt(res$ddCTRelExp[3], 100)
})

test_that("the two independent derivations of relative expression ('Johanna' and 'Henri') agree", {
  # Regression test for a sign inversion: 'Henri's ddCTPower2/ddCTRelExp previously used
  # 2^(+ddCt) instead of 2^(-ddCt), the reciprocal of the standard ddCt fold-change. Both
  # columns are meant to be two independent derivations of the same quantity.
  set.seed(7)
  Gene <- rnorm(6, mean = 25, sd = 2)
  Ref <- rnorm(6, mean = 20, sd = 1)

  res <- comparativeMethod_qPCR.analysis(Gene, Ref)

  expect_equal(res$relative_values, res$ddCTRelExp, tolerance = 1e-8)
  expect_equal(res$Power2, res$ddCTPower2, tolerance = 1e-8)
})

test_that("the calibrator sample (t0) always has a relative expression of 100", {
  Gene <- c(24, 22, 26, 21)
  Ref <- c(19, 18, 20, 17)
  res <- comparativeMethod_qPCR.analysis(Gene, Ref)

  expect_equal(res$relative_values[1], 100)
  expect_equal(res$ddCTRelExp[1], 100)
})

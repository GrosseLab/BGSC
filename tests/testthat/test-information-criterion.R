context("InformationCriterion")

test_that("BIC matches -2*logL + npar*log(k)", {
  logL <- c(-10, -8)
  npar <- c(2, 3)
  k <- c(6, 6)

  expect_equal(
    InformationCriterion(logL, npar, k, IC = "BIC"),
    -2 * logL + npar * log(k)
  )
})

test_that("AIC uses the fixed penalty of 2, independent of the number of observations k", {
  logL <- c(-10, -8)
  npar <- c(2, 3)

  expect_equal(
    InformationCriterion(logL, npar, k = c(6, 6), IC = "AIC"),
    -2 * logL + 2 * npar
  )

  # regression guard: AIC must NOT scale with k like BIC does
  expect_equal(
    InformationCriterion(logL, npar, k = c(6, 6), IC = "AIC"),
    InformationCriterion(logL, npar, k = c(1000, 1000), IC = "AIC")
  )
})

test_that("mismatched argument lengths raise an error instead of failing silently", {
  expect_error(InformationCriterion(logL = c(1, 2), npar = 1, k = c(1, 1)))
})

test_that("an unsupported IC value raises an error", {
  expect_error(InformationCriterion(logL = 1, npar = 1, k = 1, IC = "XIC"))
})

test_that("get.IC applies InformationCriterion row-wise for a matrix of log-likelihoods", {
  L <- matrix(c(-10, -12, -8, -9), nrow = 1)
  colnames(L) <- c("a", "b", "c", "d")
  npar <- c(2, 3, 3, 3)
  k <- c(6, 6, 6, 6)

  got <- get.IC(L, npar, k, IC = "BIC")
  expected <- -2 * L + matrix(npar * log(k), nrow = 1)

  expect_equal(unname(got), unname(expected))
})

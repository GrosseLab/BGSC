context("get.Posterior")

test_that("get.Posterior matches the closed-form Bayes-rule formula p(z|x) ~ exp(-BIC_z/2)*pi_z", {
  B <- matrix(c(10, 12, 15, 20,
                 5,  5,  5,  5,
                 1,  2,  3,  4),
              nrow = 3, byrow = TRUE)
  colnames(B) <- c("a", "b", "c", "d")
  Pis <- c(0.7, 0.1, 0.1, 0.1)

  got <- get.Posterior(B, Pis)

  expected <- t(apply(B, 1, function(row) {
    w <- exp(-row / 2) * Pis
    w / sum(w)
  }))
  colnames(expected) <- colnames(B)

  expect_equal(got, expected, tolerance = 1e-10)
  expect_equal(unname(rowSums(got)), rep(1, nrow(B)))
})

test_that("get.Posterior stays finite for extreme BIC values that would overflow a naive exp(-B/2)", {
  # a naive exp(-B/2) on these values overflows to Inf in double precision (B < ~-1420)
  B <- matrix(c(-3000, -2000, -1000, -500), nrow = 1)
  colnames(B) <- c("a", "b", "c", "d")

  got <- get.Posterior(B, c(0.7, 0.1, 0.1, 0.1))

  expect_true(all(is.finite(got)))
  expect_equal(sum(got), 1, tolerance = 1e-10)
  # the class with by far the smallest BIC should get essentially all posterior mass
  expect_gt(got[1, "a"], 0.999)
})

test_that("get.Posterior reduces to the prior ranking when all BIC values are equal", {
  B <- matrix(rep(0, 4), nrow = 1)
  colnames(B) <- c("a", "b", "c", "d")
  Pis <- c(0.7, 0.1, 0.1, 0.1)

  got <- get.Posterior(B, Pis)

  expect_equal(unname(got[1, ]), Pis, tolerance = 1e-10)
})

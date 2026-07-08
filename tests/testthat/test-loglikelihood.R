context("logLikelihoodOfnormData")

test_that("logLikelihoodOfnormData matches a manual Gaussian log-likelihood computation", {
  set.seed(42)
  x <- matrix(rnorm(6 * 5, mean = 8, sd = 1), nrow = 5, ncol = 6)
  rownames(x) <- paste0("g", 1:5)

  res <- logLikelihoodOfnormData(x)
  Lsets <- get.Lset()

  manualLogL <- function(d, Lset) {
    s0 <- d[Lset$s0]
    s1 <- if (!is.null(Lset$s1)) d[Lset$s1] else numeric(0)
    mu0 <- mean(s0)
    n0 <- length(s0)
    n1 <- length(s1)

    if (n1 > 0) {
      mu1 <- mean(s1)
      v <- (sum((s0 - mu0)^2) + sum((s1 - mu1)^2)) / (n0 + n1 - 2)
      resid <- c(s0 - mu0, s1 - mu1)
    } else {
      v <- sum((s0 - mu0)^2) / (n0 - 1)
      resid <- s0 - mu0
    }

    n <- n0 + n1
    -n / 2 * log(2 * pi * v) - sum(resid^2) / (2 * v)
  }

  for (g in rownames(x)) {
    for (cls in names(Lsets)) {
      expect_equal(
        unname(res$logL[g, cls]),
        manualLogL(x[g, ], Lsets[[cls]]),
        tolerance = 1e-8,
        label = paste("gene", g, "class", cls)
      )
    }
  }
})

test_that("class 'a' (no split) uses n-1 degrees of freedom for its variance estimate", {
  d <- c(1, 2, 3, 4, 5, 6)
  Lset <- get.Lset()$a
  v <- BGSC:::.Lset.get.var(d, Lset, c(s0 = mean(d), s1 = NaN), c(s0 = 6, s1 = 0))
  expect_equal(v, var(d))
})

test_that("classes b-d use a pooled n-2 variance across both groups", {
  d <- c(2, 8, 3, 9, 1, 7)
  Lset <- get.Lset()$b # s0 = c(1,3,5), s1 = c(2,4,6)
  mu0 <- mean(d[Lset$s0]); mu1 <- mean(d[Lset$s1])
  v <- BGSC:::.Lset.get.var(d, Lset, c(s0 = mu0, s1 = mu1), c(s0 = 3, s1 = 3))
  expected <- (sum((d[Lset$s0] - mu0)^2) + sum((d[Lset$s1] - mu1)^2)) / (3 + 3 - 2)
  expect_equal(v, expected)
})

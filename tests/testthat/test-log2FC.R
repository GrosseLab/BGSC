context("log2FC standard error propagation")

test_that("sqrt(SE(x)^2 + SE(y)^2) equals the Satterthwaite approximation for two independent means", {
  set.seed(3)
  x <- rnorm(20, sd = 2)
  y <- rnorm(15, sd = 3)

  combined <- sqrt(plotrix::std.error(x)^2 + plotrix::std.error(y)^2)
  satterthwaite <- sqrt(var(x) / length(x) + var(y) / length(y))

  expect_equal(combined, satterthwaite, tolerance = 1e-10)
})

test_that("comparativeMethod_qPCR.RNAi.log2FC no longer collides the C1/C0 stderr column names", {
  # Regression test for a copy-paste bug where the data.frame() calls building
  # dCT.C1C0/relVal.C1C0 used the same column name twice (intended C0 stderr column
  # was named identically to the C1 stderr column), which R silently renamed with a
  # ".1" suffix instead of erroring.
  set.seed(11)
  make_gene <- function() {
    data.frame(
      dCT = rnorm(18, mean = 5, sd = 0.5),
      relative_values = abs(rnorm(18, mean = 100, sd = 10)),
      Treatment = rep(as.character(1:6), times = 3)
    )
  }
  comp <- list(g1 = make_gene(), g2 = make_gene())

  res <- comparativeMethod_qPCR.RNAi.log2FC(comp)

  expect_true(all(c("dCT.mean.C1.stderr", "dCT.mean.C0.stderr") %in% colnames(res$dCT.C1C0)))
  expect_true(all(c("relVal.geomean.C1.stderr", "relVal.geomean.C0.stderr") %in% colnames(res$relVal.C1C0)))
  expect_false(any(grepl("\\.stderr\\.1$", colnames(res$dCT.C1C0))))
  expect_false(any(grepl("\\.stderr\\.1$", colnames(res$relVal.C1C0))))

  # the reported logFC stderr is unaffected by the naming bug: it is computed directly
  # from the correct C1/C0 stderr values before the data.frame() column-name collision
  expect_equal(
    res$dCT.C1C0$dCT.logFC.C1C0.stderr,
    sqrt(res$dCT.C1C0$dCT.mean.C1.stderr^2 + res$dCT.C1C0$dCT.mean.C0.stderr^2),
    tolerance = 1e-8
  )
})

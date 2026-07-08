context("Venn set operations")

test_that("intersect3/4/5 match base::Reduce(intersect, ...)", {
  A <- c(1, 2, 3, 4); B <- c(2, 3, 4, 5); C <- c(3, 4, 5, 6); D <- c(4, 5, 6, 7); E <- c(4, 8)

  expect_equal(intersect3(A, B, C), Reduce(intersect, list(A, B, C)))
  expect_equal(intersect4(A, B, C, D), Reduce(intersect, list(A, B, C, D)))
  expect_equal(intersect5(A, B, C, D, E), Reduce(intersect, list(A, B, C, D, E)))
})

test_that("union3/4/5 match base::Reduce(union, ...)", {
  # regression test: union5() previously called Reduce(intersect, ...) internally
  # (a copy-paste bug) and silently returned an intersection instead of a union.
  A <- c(1, 2); B <- c(2, 3); C <- c(3, 4); D <- c(4, 5); E <- c(5, 6)

  expect_equal(union3(A, B, C), Reduce(union, list(A, B, C)))
  expect_equal(union4(A, B, C, D), Reduce(union, list(A, B, C, D)))
  expect_equal(union5(A, B, C, D, E), Reduce(union, list(A, B, C, D, E)))

  # sets are disjoint enough that union and intersection give different, checkable results
  expect_equal(sort(union5(A, B, C, D, E)), 1:6)
})

test_that("dry run", {
  text <- c(
    "where the red dog lives",
    "red dog lives in the house",
    "blue cat eats mice",
    "monkeys hate cat but love trees",
    "green cat eats mice",
    "orange elephant never forgets",
    "orange elephant must forget",
    "monkeys eat banana",
    "monkeys live in trees",
    "elephant",
    "cat",
    "dog",
    "monkeys"
  )
  res <- text |>
    strsplit(" ", fixed = TRUE) |>
    gsdmm(n_iter = 100, alpha = .2, beta = .01, progress = FALSE)

  expect_true((unique(res$cluster) |> length()) < 10)
  expect_true((unique(res$cluster) |> length()) > 3)
})

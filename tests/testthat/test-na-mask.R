test_that("packed NA-mask groups match R reference across word boundaries", {
  group_na_masks <- fastphylosig:::group_na_masks_cpp
  reference_groups <- function(mask) {
    p <- ncol(mask)
    if (!p) {
      return(list(group_id = integer(), columns = list(), keep = list(),
                  count = integer(), n_keep = integer(), n_group = 0L))
    }
    labels <- vapply(seq_len(p), function(j) {
      paste(as.integer(mask[, j]), collapse = ",")
    }, character(1))
    first <- unique(labels)
    group_id <- match(labels, first)
    columns <- lapply(seq_along(first), function(i) which(group_id == i))
    keep <- lapply(seq_along(first), function(i) which(mask[, columns[[i]][1L]]))
    list(
      group_id = as.integer(group_id),
      columns = columns,
      keep = keep,
      count = lengths(columns),
      n_keep = lengths(keep),
      n_group = length(first)
    )
  }

  set.seed(801)
  for (n in c(1L, 2L, 63L, 64L, 65L, 127L, 128L, 129L, 193L)) {
    p <- sample.int(15L, 1L) + 1L
    mask <- matrix(stats::rbinom(n * p, 1, 0.55) == 1,
                   nrow = n, ncol = p)
    got <- group_na_masks(mask)
    ref <- reference_groups(mask)
    expect_equal(got$group_id, ref$group_id)
    expect_equal(got$columns, ref$columns)
    expect_equal(got$keep, ref$keep)
    expect_equal(got$count, as.integer(ref$count))
    expect_equal(got$n_keep, as.integer(ref$n_keep))
    expect_equal(got$n_group, ref$n_group)
    expect_equal(got$n_word, as.integer(ceiling(n / 64)))
    expect_length(got$key, ref$n_group)
    expect_equal(length(unique(got$key)), ref$n_group)
  }
})

test_that("all-present, all-missing, and zero-dimensional masks are coherent", {
  group_na_masks <- fastphylosig:::group_na_masks_cpp
  all_present <- matrix(TRUE, nrow = 129L, ncol = 4L)
  all_missing <- matrix(FALSE, nrow = 129L, ncol = 4L)
  for (mask in list(all_present, all_missing)) {
    got <- group_na_masks(mask)
    expect_equal(got$group_id, rep.int(1L, 4L))
    expect_equal(got$columns, list(1:4))
    expected_keep <- if (all(mask)) seq_len(nrow(mask)) else integer()
    expect_equal(got$keep, list(expected_keep))
    expect_equal(got$count, 4L)
    expect_equal(got$n_keep, length(expected_keep))
    expect_equal(got$n_group, 1L)
  }

  zero_rows <- group_na_masks(matrix(logical(), nrow = 0L, ncol = 3L))
  expect_equal(zero_rows$group_id, rep.int(1L, 3L))
  expect_equal(zero_rows$columns, list(1:3))
  expect_equal(zero_rows$keep, list(integer()))
  expect_equal(zero_rows$count, 3L)
  expect_equal(zero_rows$n_keep, 0L)
  expect_equal(zero_rows$n_group, 1L)
  expect_equal(zero_rows$n_word, 0L)

  zero_cols <- group_na_masks(matrix(logical(), nrow = 129L, ncol = 0L))
  expect_equal(zero_cols$group_id, integer())
  expect_equal(zero_cols$columns, list())
  expect_equal(zero_cols$keep, list())
  expect_equal(zero_cols$count, integer())
  expect_equal(zero_cols$n_keep, integer())
  expect_equal(zero_cols$n_group, 0L)
  expect_equal(zero_cols$n_word, 3L)
})

test_that("numeric masks are binary and NA masks are rejected", {
  group_na_masks <- fastphylosig:::group_na_masks_cpp
  numeric_mask <- matrix(c(0, 1, 1, 0, 1, 0), nrow = 3L)
  got <- group_na_masks(numeric_mask)
  expect_equal(got$group_id, c(1L, 2L))
  expect_equal(got$keep, list(c(2L, 3L), 2L))
  expect_error(group_na_masks(matrix(c(TRUE, NA), nrow = 2L)),
               "contains NA")
  expect_error(group_na_masks(matrix(c(0, NA_real_), nrow = 2L)),
               "contains NA")
  expect_error(group_na_masks(matrix(c(0, 2), nrow = 2L)),
               "only 0 and 1")
})

test_that("different packed words never merge distinct patterns", {
  group_na_masks <- fastphylosig:::group_na_masks_cpp
  n <- 193L
  mask <- matrix(FALSE, nrow = n, ncol = 5L)
  mask[c(1L), 2L] <- TRUE
  mask[c(65L), 3L] <- TRUE
  mask[c(129L), 4L] <- TRUE
  mask[c(1L, 65L, 129L), 5L] <- TRUE
  got <- group_na_masks(mask)
  expect_equal(got$group_id, 1:5)
  expect_equal(got$columns, lapply(1:5, identity))
  expect_equal(got$keep, list(integer(), 1L, 65L, 129L, c(1L, 65L, 129L)))
  expect_equal(got$count, rep.int(1L, 5L))
})

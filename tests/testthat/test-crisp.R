set.seed(123)
# data
xu <- stats::rnorm(50)
xm <- matrix(stats::rnorm(100), ncol = 2)
index <- seq(Sys.Date(), Sys.Date() + 49, "day")
cond <- xu < 0

# Univariate
xu_vec <- xu
xu_mat <- as.matrix(xu)
xu_ts <- stats::as.ts(xu)
xu_xts <- xts::xts(xu, index)
xu_df <- as.data.frame(xu)
xu_tbl <- tibble::tibble(index = index, x = xu)

# Multivariate
xm_vec <- xm
xm_mat <- as.matrix(xm)
xm_ts <- stats::as.ts(xm)
xm_xts <- xts::xts(xm, index)
xm_df <- as.data.frame(xm)
xm_tbl <- tibble::tibble(index = index, x = xm)

# condition must be specified ---------------------------------------------

test_that("condition must specified", {
  expect_error(crisp(xu))
})

test_that("error if condition is not logical", {
  expect_error(crisp(xu, as.data.frame(cond)))
  expect_error(crisp(xu, as.matrix(cond)))
})

test_that("error if sizes differ", {
  expect_error(crisp(c(xu, rnorm(2)), cond))
})


# works on different classes ----------------------------------------------

# doubles
crisp_numeric_u <- crisp(xu, cond)
test_that("works on univariate doubles", {
  # type
  expect_type(crisp_numeric_u, "double")
  expect_s3_class(crisp_numeric_u, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_numeric_u), vctrs::vec_size(xu))

})


# matrices
crisp_matu <- crisp(xu_mat, cond)
test_that("works on univariate matrices", {
  # type
  expect_type(crisp_matu, "double")
  expect_s3_class(crisp_matu, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_matu), vctrs::vec_size(xu_mat))
})

crisp_matm <- crisp(xm_mat, cond)
test_that("works on multivariate matrices", {
  # type
  expect_type(crisp_matm, "double")
  expect_s3_class(crisp_matm, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_matm), vctrs::vec_size(xm_mat))
})

# ts
crisp_tsu <- crisp(xu_ts, cond)
test_that("works on univariate ts", {
  # type
  expect_type(crisp_tsu, "double")
  expect_s3_class(crisp_tsu, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_tsu), vctrs::vec_size(xu_ts))
})

crisp_tsm <- crisp(xm_ts, cond)
test_that("works on multivariate ts", {
  # type
  expect_type(crisp_tsm, "double")
  expect_s3_class(crisp_tsm, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_tsm), vctrs::vec_size(xm_ts))
})

# xts
crisp_xtsu <- crisp(xu_xts, cond)
test_that("works on univariate xts", {
  # type
  expect_type(crisp_xtsu, "double")
  expect_s3_class(crisp_xtsu, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_xtsu), vctrs::vec_size(xu_xts))
})

crisp_xtsm <- crisp(xm_xts, cond)
test_that("works on multivariate xts", {
  # type
  expect_type(crisp_xtsm, "double")
  expect_s3_class(crisp_xtsm, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_xtsm), vctrs::vec_size(xm_xts))
})

# data.frame
crisp_dfu <- crisp(xu_df, cond)
test_that("works on univariate data.frames", {
  # type
  expect_type(crisp_dfu, "double")
  expect_s3_class(crisp_dfu, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_dfu), vctrs::vec_size(xu_df))
})

crisp_dfm <- crisp(xm_df, cond)
test_that("works on multivariate data.frames", {
  # type
  expect_type(crisp_dfm, "double")
  expect_s3_class(crisp_dfm, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_dfm), vctrs::vec_size(xm_df))
})

# tbl
crisp_tblu <- crisp(xu_tbl, cond)
test_that("works on univariate tibbles", {
  # type
  expect_type(crisp_tblu, "double")
  expect_s3_class(crisp_tblu, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_tblu), vctrs::vec_size(xu_tbl))
})

crisp_tblm <- crisp(xm_tbl, cond)
test_that("works on multivariate tibbles", {
  # type
  expect_type(crisp_tblm, "double")
  expect_s3_class(crisp_tblm, "ffp")
  # size
  expect_equal(vctrs::vec_size(crisp_tblm), vctrs::vec_size(xm_tbl))
})

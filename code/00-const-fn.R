## 00-const-fn.R -----------------------------------------------------
##
## Define constants and helper functions.


### Function definitions ---------------------------------------------

#' Produce a data set for use as an input for effect plots for the
#' main SN analysis.
#'
#' @param sizes Character vector of breed sizes. Default `NULL` means
#'   all sizes.
#' @param sexes Character vector of sexes. Default `NULL` means all
#'   sexes.
#' @param ages Numeric vector of ages. Default `NULL` means all ages
#'   from 0.5 years to `max(ageYearsX)` years in increments of 0.5
#'   years.
#' @param weight_pctls Numeric vector of weight percentiles. Default
#'   `NULL` means 25, 50, and 75, i.e., 1st quartile, median, and
#'   third quartile.
#' @return A data.frame with a row for each combination of the input
#'   values.
#' @examples
#' define_sn_reference_points()
#' define_sn_reference_points(
#'   ages = seq(0.5, 6.5, by = 0.5), weight_pctls = 50)
define_sn_reference_points <- function(
  sizes = NULL, sexes = NULL, ages = NULL, weight_pctls = NULL
) {
  # Define a helper function to get the actual values of weight for
  # the requested percentile.
  get_wt_qntl <- function(x) {
    # For integer values of age, use records of this size and sex and
    # where rounded age equals the reference value. For noninteger
    # values of age, use records of this size and sex and where
    # rounded age equals the floor or ceiling of the reference value;
    # then average the weights at the quantiles of the two ages.
    if (as.double(x[["age"]]) %% 1 == 0) {
      dat[
        ageYearsR == as.double(x[["age"]]) & size == x[["size"]]
        & sex == x[["sex"]],
        quantile(weight, as.double(x[["wt_pctl"]]) / 100)]
    } else {
      mean(
        c(dat[
            ageYearsR == floor(as.double(x[["age"]]))
            & size == x[["size"]] & sex == x[["sex"]],
            quantile(weight, as.double(x[["wt_pctl"]]) / 100)],
          dat[
            ageYearsR == ceiling(as.double(x[["age"]]))
            & size == x[["size"]] & sex == x[["sex"]],
            quantile(weight, as.double(x[["wt_pctl"]]) / 100)]))
    }
  }

  if (!is.null(sizes)) {
    size_ref_pts <- sizes
  } else {
    size_ref_pts <- levels(dat$size)
  }
  if (!is.null(sexes)) {
    sex_ref_pts <- sexes
  } else {
    sex_ref_pts <- levels(dat$sex)
  }
  if (!is.null(ages)) {
    age_ref_pts <- ages
  } else {
    age_ref_pts <- seq(0.5, max(dat$ageYearsX), by = 0.5)
  }
  if (!is.null(weight_pctls)) {
    wt_pctl <- weight_pctls
  } else {
    wt_pctl <- c(25, 50, 75)
  }

  ref_ds <- expand.grid(
    size = size_ref_pts, sex = sex_ref_pts, age = age_ref_pts,
    wt_pctl = wt_pctl)
  ref_ds$weight <- apply(ref_ds, 1, get_wt_qntl)
  ref_ds
}


#' Produce a data set for use as an input for effect plots for the age
#' effect among SN analysis.
#'
#' @param sizes Character vector of breed sizes. Default `NULL` means
#'   all sizes.
#' @param sexes Character vector of sexes. Default `NULL` means all
#'   sexes.
#' @param ages Numeric vector of ages. Default `NULL` means all ages
#'   from 0.5 years to `max(ageYearsX)` years in increments of 0.5
#'   years.
#' @return A data.frame with a row for each combination of the input
#'   values.
#' @examples
#' define_age_reference_points()
#' define_age_reference_points(ages = seq(0.5, 6.5, by = 0.5))
define_age_reference_points <- function(
  sizes = NULL, sexes = NULL, ages = NULL
) {
  if (!is.null(sizes)) {
    size_ref_pts <- sizes
  } else {
    size_ref_pts <- levels(dat$size)
  }
  if (!is.null(sexes)) {
    sex_ref_pts <- sexes
  } else {
    sex_ref_pts <- levels(dat$sex)
  }
  if (!is.null(ages)) {
    age_ref_pts <- ages
  } else {
    age_ref_pts <- seq(0.5, max(dat$ageYearsX), by = 0.5)
  }

  ref_ds <- expand.grid(
    comparator_age = age_ref_pts, reference_age = age_ref_pts,
    sex = sex_ref_pts, size = size_ref_pts)
  ref_ds[, c("size", "sex", "reference_age", "comparator_age")]
}


evaluate_sn_reference_points <- function(reference_points, model) {
  eval_at <- function(x) {
    tbl <- summary(
      model,
      sex = x[["sex"]],
      ageYearsX = as.numeric(x[["age"]]),
      weight = as.numeric(x[["weight"]])
    )
    sn_row <- which(rownames(tbl) == "sn - Spayed/neutered:Intact")
    hr <- tbl[sn_row+1, "Effect"]
    lo <- tbl[sn_row+1, "Lower 0.95"]
    hi <- tbl[sn_row+1, "Upper 0.95"]
    c(hr = hr, lo = lo, hi = hi)
  }

  res <- apply(reference_points, 1, eval_at)
  t(res)
}


evaluate_age_among_sn_reference_points <- function(
  reference_points, model
) {
  eval_at <- function(x) {
    tbl <- summary(
      model,
      sn = "Spayed/neutered",
      sex = x[["sex"]],
      ageYearsX = c(
        as.numeric(x[["reference_age"]]),
        as.numeric(x[["comparator_age"]]))
    )
    age_row <- which(rownames(tbl) == "ageYearsX")
    hr <- tbl[age_row+1, "Effect"]
    lo <- tbl[age_row+1, "Lower 0.95"]
    hi <- tbl[age_row+1, "Upper 0.95"]
    c(hr = hr, lo = lo, hi = hi)
  }

  res <- apply(reference_points, 1, eval_at)
  t(res)
}

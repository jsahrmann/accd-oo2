## 03-oo-coxph.R -----------------------------------------------------
##
## Fit a Cox proportional hazards model as part of the primary
## analysis for the ACC&D overweight/obesity analysis.


### Setup ------------------------------------------------------------

library(data.table)
library(Hmisc)
library(magrittr)
library(rms)
library(survival)


### Input data -------------------------------------------------------

source("./00-const-fn.R")
source("./02-read-cohort.R")


### Data management --------------------------------------------------

sort(
  table(dat[size == "Toy and Small", breed]), decreasing = TRUE)[1:10]
sort(table(dat[size == "Large", breed]), decreasing = TRUE)[1:5]

toysmall_dat <- dat[size == "Toy and Small"]
york_dat <- toysmall_dat[breed == "Yorkshire Terrier"]
sansyork_dat <- toysmall_dat[breed != "Yorkshire Terrier"]

Hmisc::units(toysmall_dat$t2e_oo, "day")
dd <- datadist(toysmall_dat)
options(datadist = "dd")

surv_oo <- toysmall_dat %$% Surv(oo_t2e, oo_event)
toysmall_mod <- cph(
  surv_oo ~
    sn + rcs(ageYearsX, 3) + sex + mixed_breed +
    wellness_plan + rcs(weight, 3) + rcs(visitsPerYear, 3) +
    sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
    sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
  data = toysmall_dat, x = TRUE, y = TRUE, surv = TRUE)


large_dat <- dat[size == "Large"]
pug_dat <- large_dat[breed == "Pug"]
sanspug_dat <- large_dat[breed != "Pug"]
Hmisc::units(sanspug_dat$oo_t2e, "day")
dd <- datadist(sanspug_dat)
options(datadist = "dd")
surv_oo <- sanspug_dat %$% Surv(oo_t2e, oo_event)
sanspug_oo_mod <- cph(
  surv_oo ~
    sn + rcs(ageYearsX, 3) + sex + mixed_breed +
    wellness_plan + rcs(weight, 3) + rcs(visitsPerYear, 3) +
    sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
    sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
  data = sanspug_dat, x = TRUE, y = TRUE, surv = TRUE)
Hmisc::units(sanspug_dat$obese_t2e, "day")
dd <- datadist(sanspug_dat)
options(datadist = "dd")
surv_ob <- sanspug_dat %$% Surv(obese_t2e, oo_event)
sanspug_ob_mod <- cph(
  surv_ob ~
    sn + rcs(ageYearsX, 5) + sex + mixed_breed +
    wellness_plan + rcs(weight, 3) + rcs(visitsPerYear, 3) +
    sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
    sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
  data = sanspug_dat, x = TRUE, y = TRUE, surv = TRUE)

oo2_make_data_sets <- function(.data, .size, .breed) {
  size_dat <- .data[size == .size]
  breed_dat <- size_dat[breed == .breed]
  sans_breed_dat <- size_dat[breed != .breed]
  list(
    size_dat = size_dat,
    breed_dat = breed_dat,
    sans_breed_dat = sans_breed_dat)
}

oo2_fit_models <- function(ds_list) {
  # Define data distributions for each data set.
  size_datadist <<- rms::datadist(ds_list$size_dat)
  breed_datadist <- rms::datadist(ds_list$breed_dat)
  sans_breed_datadist <- rms::datadist(ds_list$sans_breed_dat)
  # Fit the O/O model for each data set.
  options(datadist = "size_datadist")
  size_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$size_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  options(datadist = "breed_datadist")
  breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  options(datadist = "sans_breed_datadist")
  sans_breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$sans_breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  # Fit the obese model for each data set.
  options(datadist = "size_datadist")
  size_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$size_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  options(datadist = "breed_datadist")
  breed_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  options(datadist = "sans_breed_datadist")
  sans_breed_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$sans_breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  # Collect results into a list.
  list(
    size = list(
      datadist = size_datadist,
      oo_mod = size_oo_mod,
      ob_mod = size_ob_mod
    ),
    breed = list(
      datadist = breed_datadist,
      oo_mod = breed_oo_mod,
      ob_mod = breed_ob_mod
    ),
    sans_breed = list(
      datadist = sans_breed_datadist,
      oo_mod = sans_breed_oo_mod,
      ob_mod = sans_breed_ob_mod
    )
  )
}

oo2_fit_models <- function(ds_list) {
  # Define data distributions for each data set.
  size_datadist <<- rms::datadist(ds_list$size_dat)
  breed_datadist <<- rms::datadist(ds_list$breed_dat)
  sans_breed_datadist <<- rms::datadist(ds_list$sans_breed_dat)
  # Fit the O/O model for each data set.
  size_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$size_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  sans_breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$sans_breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  # Fit the obese model for each data set.
  size_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$size_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  breed_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  sans_breed_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 3) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 3) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 3) + sex : rcs(weight, 3),
    data = ds_list$sans_breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  # Collect results into a list.
  list(
    size = list(
      oo_mod = size_oo_mod,
      ob_mod = size_ob_mod
    ),
    breed = list(
      oo_mod = breed_oo_mod,
      ob_mod = breed_ob_mod
    ),
    sans_breed = list(
      oo_mod = sans_breed_oo_mod,
      ob_mod = sans_breed_ob_mod
    )
  )
}

pug_ds <- oo2_make_data_sets(dat, "Toy and Small", "Pug")
pug_models <- oo2_fit_models(pug_ds)

sn_ref_pts <- data.table::as.data.table(
  define_sn_reference_points(
    sizes = "Toy and Small",
    ages = c(0.25, seq(0.5, 6, by = 0.5))
  )
)

options(datadist = "size_datadist")
size_oo_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$size$oo_mod)
size_ob_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$size$ob_mod)
options(datadist = "breed_datadist")
breed_oo_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$breed$oo_mod)
breed_ob_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$breed$ob_mod)
options(datadist = "sans_breed_datadist")
sans_breed_oo_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$sans_breed$oo_mod)
sans_breed_ob_sn_est <- evaluate_sn_reference_points(
  sn_ref_pts, model = pug_models$sans_breed$ob_mod)


### Modeling ---------------------------------------------------------

# Set options for Hmisc/rms analysis functions.
Hmisc::units(dat$t2e_oo, "day")
dd <- datadist(dat)
options(datadist = "dd")

# Create the survival outcome object.
surv_oo <- dat %$% Surv(oo_t2e, oo_event)

# Fit the primary model.
model1 <- cph(
  surv_oo ~
    sn + rcs(ageYearsX, 5) + size + sex + mixed_breed +
    wellness_plan + rcs(weight, 3) + rcs(visitsPerYear, 3) +
    sn : rcs(ageYearsX, 5) + sn : size + sn : sex + sn : rcs(weight, 3) +
    size : rcs(ageYearsX, 5) + size : rcs(weight, 3) + size : sex +
    sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
  data = dat, x = TRUE, y = TRUE, surv = TRUE)

anova(model1)
summary(model1)

# Get hazard ratio for one-unit increase in average number of visits
# per year.
summary(model1, visitsPerYear = 1:2)


# SN effect -----------------------

sn_ref_pts <- as.data.table(
  define_sn_reference_points(ages = c(0.25, seq(0.5, 6, by = 0.5))))

sn_est <- evaluate_sn_reference_points(sn_ref_pts, model1)

sn_results <- cbind(sn_ref_pts, sn_est)


# Age effect among SN -------------

age_among_sn_ref_pts <- data.table::as.data.table(
  define_age_reference_points(ages = c(0.25, seq(0.5, 6, by = 0.5))))

age_among_sn_est <- evaluate_age_among_sn_reference_points(
  age_among_sn_ref_pts, model1)

age_among_sn_results <- cbind(age_among_sn_ref_pts, age_among_sn_est)


# Exports ------------------------------------------------------------

save(
  sn_results, age_among_sn_results,
  file = "../data/oo-results.Rdata"
)


# SN effect -----------------------

# Smallest difference in risk by SN
sn_results[wt_pctl == 50][which.min(hr)]
# Largest difference in risk by SN
sn_results[wt_pctl == 50][which.max(hr)]

round(sn_results[wt_pctl == 50][which.max(hr)]$hr, digits = 2)
round(sn_results[wt_pctl == 50][which.max(hr)]$hi, digits = 2)


# Age effect among SN -------------

# Smallest difference in risk by SN
age_among_sn_results[reference_age == 1.0][which.min(hr)]

# Largest difference in risk by SN
age_among_sn_results[wt_pctl == 50][which.max(hr)]

round(age_among_sn_results[wt_pctl == 50][which.max(hr)]$hr, digits = 2)
round(age_among_sn_results[wt_pctl == 50][which.max(hr)]$hi, digits = 2)

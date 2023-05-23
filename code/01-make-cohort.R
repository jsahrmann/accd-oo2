## 01-make-cohort.R --------------------------------------------------
##
## Data management for the ACC&D analysis of the effect of
## sterilization on risk of overweight/obese status.


### Setup ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggconsort)
library(lubridate)
library(magrittr)
library(readxl)


### Constants --------------------------------------------------------

data_dir <- "~/Dropbox/Banfield Dog Data/Data Files from Banfield"

new_col_names <- c(
  "id", "sex", "birth_date", "breed", "mixed_breed", "neuter_date",
  "first_visit_date", "first_overweight_date", "state", "visit_date",
  "intact", "visit_reason", "wellness_plan", "bcs", "weight",
  "dx_ccl", "dx_hypothyroidism", "dx_hyperthyroidism")


### Input data -------------------------------------------------------

# Read the full Banfield data set.
visits_all <- data.table::fread(
  file = paste(data_dir, "R2020_ACCD_DATA.txt", sep = "/"),
  na.strings = "null"
)
data.table::setnames(visits_all, new_col_names)

# Recast columns with dates.
date_cols <- grep("date", colnames(visits_all), value = TRUE)
visits_all[, (date_cols) := lapply(.SD, mdy), .SDcols = date_cols]

# Read the breed sizes data set.

breed_sizes <- readxl::read_excel(
  "~/home/accd-oo/data/Banfield data breeds - final recommendations.xlsx",
  range = "A1:I506"
) %>%
  data.table::as.data.table()

data.table::setnames(
  breed_sizes,
  old = c(
    "Banfield breed label",
    "Recommended size category (DAP buckets)"
  ),
  new = c("breed", "size")
)

breed_sizes <- breed_sizes[, .(breed, size)]


### Data management --------------------------------------------------

# Define an index date for each dog as `neuter_date` for dogs who were
# spayed/neutered in 2014 or the earliest visit date in 2014 for those
# who were not. Only consider dogs who were intact as of December 31,
# 2013. (Dogs who were overweight/obese in 2013 have already been
# excluded by Banfield.)
visits_2014ButNotIfSNIn2013 <- visits_all[
  lubridate::year(visit_date) == 2014
  & (is.na(neuter_date) | lubridate::year(neuter_date) > 2013)
][,
  `:=`(
    sn = data.table::fifelse(
      !is.na(neuter_date) & lubridate::year(neuter_date) == 2014,
      1, 0),
    index_date = data.table::fifelse(
      !is.na(neuter_date) & lubridate::year(neuter_date) == 2014,
      neuter_date, min(visit_date))
  ),
  by = id
]

# Make a data set of unique dogs in this sample for later merging.
dogs_a14 <- visits_2014ButNotIfSNIn2013[,
  lapply(.SD, first), by = id, .SDcols = c("sn", "index_date")
]

# Select all visits for dogs in this sample.
visits_a14_all <- merge(dogs_a14, visits_all, by = "id")

# Add breed size category to the 'all visits' data set.
visits_a14_all <- merge(
  visits_a14_all, breed_sizes, by = "breed", all.x = TRUE)


# Exclusions ---------------------------------------------------------


# At or before index --------------
#
# - Overweight or obese any time at or before index
# - Hypothyroidism diagnosis any time at or before index
# - Hyperthyroidism diagnosis any time at or before index

# Select visits at or before index.
visits_a14_atIndexOrBefore <- visits_a14_all[visit_date <= index_date]

visits_a14_atIndexOrBefore[,
  `:=`(
    first_visit_date_ever = first_visit_date,
    first_visit_date = min(visit_date)
  ),
  by = id
][,
  dt2013VisitAndIndex := as.integer(index_date - first_visit_date + 1)
][,
  nVisits := dplyr::n_distinct(visit_date), by = id
][,
  visitsPerYear := nVisits / dt2013VisitAndIndex * 365.25
]

# Filter to visits with a BCS of 4 or 5.
visits_a14OO_atIndexOrBefore <- visits_a14_atIndexOrBefore[
  bcs %in% 4:5]

# Collapse to one record per dog.
dogs_a14ExclOO_atIndexOrBefore <- unique(
  visits_a14OO_atIndexOrBefore[, .(id, sn, index_date)])

# For each dog, check for a history of conditions that could be
# related to weight.
dogs_a14_atIndexOrBefore <- visits_a14_atIndexOrBefore[,
  .(dx_hypothy_preEver = max(dx_hypothyroidism),
    dx_hyperthy_preEver = max(dx_hyperthyroidism)),
  by = .(id, sn)]

# Make a data set of IDs of dogs with a weight-related exclusion
# diagnosis at or before index.
dogs_a14ExclDx_atIndexOrBefore <- dogs_a14_atIndexOrBefore[
    dx_hypothy_preEver == 1 | dx_hyperthy_preEver == 1, .(id, sn)]


# At index ------------------------

# - neuter date but no corresponding visit
# - neuter date missing but `intact == 0`
# - unrealistically high value for weight
# - BCS == 1
# - unresolved breed size
# - S/N at less than 90 days of age

# Subset to visits at index.
visits_a14_atIndex <- visits_a14_all[visit_date == index_date]

# Count the number of visit records on the index date for each dog.
dogs_a14VisitCounts_atIndex <- visits_a14_atIndex[,
  .(visits = .N), by = id]

# To simplify, we'll take the last record per dog. (Some hacking
# suggests that the difference between taking the first or last is
# minimal.)
dogs_a14_atIndex <- visits_a14_atIndex[, lapply(.SD, last), by = id]

# A small number of dogs don't appear in the data set at index. These
# dogs all have a neuter date (which is defined in the data dictionary
# as the spay/neuter date *at Banfield*) but no corresponding visit
# record.
dogs_a14ExclSNButNoVisit_atIndex <- dogs_a14[
  !dogs_a14_atIndex, on = "id"]

# Select dogs with missing BCS or BCS of 1 or with an abnormal weight.
dogs_a14ExclBCSMissing_atIndex <- dogs_a14_atIndex[
  is.na(bcs), .(id, sn)]
dogs_a14ExclBCS1_atIndex <- dogs_a14_atIndex[
  bcs == 1, .(id, sn)]
dogs_a14ExclWt_atIndex <- dogs_a14_atIndex[
  weight >= 250, .(id, sn)]

# Select dogs with missing neuter date but `intact == 0`, indicating
# spayed/neutered at the start of the index visit. (Perhaps these dogs
# were sterilized elsewhere? Curiously, `sex` does not suggest
# spayed/neutered.)
dogs_a14ExclNoNeuterDateButNotIntact_atIndex <- dogs_a14_atIndex[
  is.na(neuter_date) & intact == 0, .(id, sn)]

# Select dogs whose breed could not be categorized by size (or that
# are otherwise unusual).
dogs_a14ExclBreed_atIndex <- dogs_a14_atIndex[is.na(size), .(id, sn)]

# Select spayed/neutered dogs who were sterilized dogs at less than 90
# days of age.
dogs_a14ExclEarlySN_atIndex <- dogs_a14_atIndex[,
  age := as.integer(index_date - birth_date)
][
  sn == 1 & age < 90, .(id, sn)
]


# After index ---------------------

# Select dogs with no visits after index for 'exclusion', as they won't
# contribute any information to the survival analysis.
dogs_a14WithFoo <- unique(
  visits_a14_all[visit_date > index_date, .(id, sn)])
dogs_a14WithoutFoo <- dogs_a14[!dogs_a14WithFoo, .(id, sn), on = "id"]


# Final index data set -----------------------------------------------

# Apply exclusions in the order recommended by Jan in her 20220522
# revision of the Results.
#
# "Keep: No follow-up visits"
dogs_i1 <- dogs_a14[!dogs_a14WithoutFoo, .(id, sn), on = "id"]
dogs_x1 <- dogs_a14[!dogs_i1, .(id, sn), on = "id"]
# "Keep: No visit record on S/N date"
dogs_i2 <- dogs_i1[
  !dogs_a14ExclSNButNoVisit_atIndex, .(id, sn), on = "id"]
dogs_x2 <- dogs_i1[!dogs_i2, .(id, sn), on = "id"]
# "Combine: S/N before 90 days with neutered elsewhere"
dogs_i3 <- dogs_i2[
  !dogs_a14ExclEarlySN_atIndex, .(id, sn), on = "id"
][
  !dogs_a14ExclNoNeuterDateButNotIntact_atIndex, .(id, sn), on = "id"
]
dogs_x3 <- dogs_i2[!dogs_i3, .(id, sn), on = "id"]
# "Keep: BCS missing at index date"
dogs_i4 <- dogs_i3[
  !dogs_a14ExclBCSMissing_atIndex, .(id, sn), on = "id"]
dogs_x4 <- dogs_i3[!dogs_i4, .(id, sn), on = "id"]
# "Combine: BCS < 3 or >3 at or before index date or weight 250
# lbs. or more"
dogs_i5 <- dogs_i4[
  !dogs_a14ExclBCS1_atIndex, .(id, sn), on = "id"
][
  !dogs_a14ExclOO_atIndexOrBefore, .(id, sn), on = "id"
][
  !dogs_a14ExclWt_atIndex, .(id, sn), on = "id"
]
dogs_x5 <- dogs_i4[!dogs_i5, .(id, sn), on = "id"]
# "Keep: Hyperthyroid or hypothyroid at or before index date"
dogs_i6 <- dogs_i5[
  !dogs_a14ExclDx_atIndexOrBefore, .(id, sn), on = "id"]
dogs_x6 <- dogs_i5[!dogs_i6, .(id, sn), on = "id"]
# "Keep: Size category could not be assigned"
dogs_i7 <- dogs_i6[!dogs_a14ExclBreed_atIndex, .(id, sn), on = "id"]
dogs_x7 <- dogs_i6[!dogs_i7, .(id, sn), on = "id"]

# Print the sample size changes for a flow chart table.
purrr::map(
  list(
    dogs_a14, dogs_x1, dogs_i1, dogs_x2, dogs_i2, dogs_x3, dogs_i3,
    dogs_x4, dogs_i4, dogs_x5, dogs_i5, dogs_x6, dogs_i6, dogs_x7,
    dogs_i7
  ),
  ~ table(.x[["sn"]])
)

# Select the final sample for inclusion.
dogs_final <- dogs_a14[dogs_i7, on = c("id", "sn")]

# Add characteristics from index visits.
dogs_final <- dogs_a14_atIndex[dogs_final[, .(id)], on = "id"]

# Compute the rate of visits per year based on the pre-index and index
# periods.
visits_final <- visits_a14_all[dogs_final[, .(id)], on = "id"]
visits_final_atIndexOrBefore <- visits_final[visit_date <= index_date]

visits_final_atIndexOrBefore[,
  `:=`(
    first_visit_date_ever = first_visit_date,
    first_visit_date = min(visit_date)
  ),
  by = id
][,
  dt2013VisitAndIndex := as.integer(index_date - first_visit_date + 1)
][,
  nVisits := n_distinct(visit_date), by = id
][,
  visitsPerYear := nVisits / dt2013VisitAndIndex * 365.25
]

dogs_finalVisitsPerYear <- visits_final_atIndexOrBefore[,
  lapply(.SD, data.table::first), by = id, .SDcols = "visitsPerYear"
]

# Add `visitsPerYear` to the final inclusion data set.
dogs_final <- dogs_final[dogs_finalVisitsPerYear, on = "id"]

# Define new variables at index and discard unneeded ones.
dogs_final[,
  `:=`(
    ageDays = age,
    sex = fcase(
      sex == "Neutered Male", "Male",
      sex == "Spayed Female", "Female",
      sex == "Male", "Male",
      sex == "Female", "Female"),
    lenpre = as.integer(index_date - first_visit_date)
  )
][,
  c(
    "birth_date", "neuter_date", "first_visit_date",
    "first_overweight_date", "visit_date", "intact", "bcs", "dx_ccl",
    "dx_hypothyroidism", "dx_hyperthyroidism", "age")
  := NULL
]

# Per Valerie 20230510-2022, combine sub-breeds for chihuahuas,
# dachshunds, and poodles.
dogs_final[,
  `:=`(
    breed = data.table::fifelse(
      breed %in% c(
        "Poodle", "Toy Poodle", "Miniature Poodle", "Teacup Poodle",
        "French Poodle", "Medium Poodle"
      ),
      "Poodle",
      breed
    )
  )
][,
  `:=`(
    breed = data.table::fifelse(
      breed %in% c(
        "Chihuahua", "Chihuahua (Smooth Coat)",
        "Chihuahua (Long Coat)"
      ),
      "Chihuahua",
      breed
    )
  )
][,
  `:=`(
    breed = data.table::fifelse(
      breed %in% c(
        "Dachshund", "Long Haired Dachshund",
        "Miniature Long Haired Dachshund", "Short Haired Dachshund",
        "Miniature Wire Haired Dachshund", "Smooth Haired Dachshund",
        "Wirehair Dachshund", "Rabbit Dachshund"
      ),
      "Dachshund",
      breed
    )
  )
]

readr::write_rds(
  dogs_final,
  paste0(
    "~/Dropbox/Banfield Dog Data/R Data Files (for analysis)/",
    "dogs_final_20230522.rds"
  ),
  compress = "gz"
)


### Outcomes and censoring -------------------------------------------

# Make a data set of post-index visit records.
visits_final_postIndex <-
  visits_a14_all[dogs_final[, .(id)], on = "id"]

# Define dates for outcome and censoring events. (Note we're using the
# Date type-specific form of the missing value provided by the
# lubridate package, otherwise we get a type mismatch error.)
visits_final_postIndex[,
  `:=`(
    neuter_date_after2014 = fifelse(
      !is.na(neuter_date) & year(neuter_date) > 2014,
      neuter_date, NA_Date_),
    dx_hypothy_date = fifelse(
      dx_hypothyroidism == 1, visit_date, NA_Date_),
    dx_hyprthy_date = fifelse(
      dx_hyperthyroidism == 1, visit_date, NA_Date_),
    oo_date = fifelse(
      bcs == 4 | bcs == 5, visit_date, NA_Date_),
    obese_date = fifelse(
      bcs == 5, visit_date, NA_Date_))]

# Choose the earliest date of each outcome or censoring event per dog.
dogs_outcomeAndCensoringDates <- visits_final_postIndex[,
  .(last_visit_date = max(visit_date),
    neuter_date_after2014 = min(neuter_date_after2014, na.rm = TRUE),
    dx_hypothy_date_earliest = min(dx_hypothy_date, na.rm = TRUE),
    dx_hyprthy_date_earliest = min(dx_hyprthy_date, na.rm = TRUE),
    oo_date_earliest = min(oo_date, na.rm = TRUE),
    obese_date_earliest = min(obese_date, na.rm = TRUE)),
  by = .(id, index_date)]

# For dogs who never experience a certain type of outcome or censoring
# event (e.g., dogs never diagnosed with hypothyroidism), the
# aggregation functions return `Inf` rather than `NA`, so we need to
# manually set these to missing.
dogs_outcomeAndCensoringDates[,
  `:=`(
    dx_hypothy_date_earliest = fifelse(
      is.infinite(dx_hypothy_date_earliest),
      NA_Date_,
      dx_hypothy_date_earliest),
    dx_hyprthy_date_earliest = fifelse(
      is.infinite(dx_hyprthy_date_earliest),
      NA_Date_,
      dx_hyprthy_date_earliest),
    oo_date_earliest = fifelse(
      is.infinite(oo_date_earliest),
      NA_Date_,
      oo_date_earliest),
    obese_date_earliest = fifelse(
      is.infinite(obese_date_earliest),
      NA_Date_,
      obese_date_earliest))
][,  # Choose the earliest event per outcome.
  `:=`(
    oo_event_date = pmin(
      last_visit_date, neuter_date_after2014,
      dx_hypothy_date_earliest, dx_hyprthy_date_earliest,
      oo_date_earliest,
      na.rm = TRUE),
    obese_event_date = pmin(
      last_visit_date, neuter_date_after2014,
      dx_hypothy_date_earliest, dx_hyprthy_date_earliest,
      obese_date_earliest,
      na.rm = TRUE))
][,  # Define event indicator and time-to-event variables.
  `:=`(
    oo_event = fcase(
      is.na(oo_event_date), 0,
      oo_date_earliest == oo_event_date, 1,
      default = 0),
    obese_event = fcase(
      is.na(obese_event_date), 0,
      obese_date_earliest == obese_event_date, 1,
      default = 0),
    oo_t2e = as.integer(oo_event_date - index_date),
    obese_t2e = as.integer(obese_event_date - index_date))
]

readr::write_rds(
  dogs_outcomeAndCensoringDates,
  paste0(
    "~/Dropbox/Banfield Dog Data/R Data Files (for analysis)/",
    "dogs_outcomeAndCensoringDates_20230522.rds"
  ),
  compress = "gz"
)

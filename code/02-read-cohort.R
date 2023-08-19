## 02-read-cohort.R --------------------------------------------------
##
## Read the analytic data set produced by 01-make-cohort.R.


### Setup ------------------------------------------------------------

library(data.table)


### Input data -------------------------------------------------------

# Read and join the covariates and outcomes data sets.
basedat <- readr::read_rds(
  paste0(
    "~/Dropbox/Banfield Dog Data/R Data Files (for analysis)/",
    "dogs_final_20230522.rds")
)
survdat <- readr::read_rds(
  paste0(
    "~/Dropbox/Banfield Dog Data/R Data Files (for analysis)/",
    "dogs_outcomeAndCensoringDates_20230819.rds")
)
dat <- basedat[survdat, on = c("id", "index_date")]
rm(basedat, survdat)


### Data management --------------------------------------------------

# Define additional variables for analysis, and declare factor
# variables.
dat[,
  `:=`(
    ageYearsX = ageDays / 365.25,
    ageYearsT = floor(ageDays / 365.25),
    ageYearsR = round(ageDays / 365.25),
    sex = factor(sex, levels = c("Female", "Male")),
    mixed_breed = factor(mixed_breed, levels = c("N", "Y")),
    size = factor(
      size, levels = c(
        "Toy and Small", "Medium", "Standard", "Large", "Giant")),
    sn = factor(
      sn, levels = 0:1, labels = c("Intact", "Spayed/neutered"))    
  )
]

# Cap follow-up at five years.
dat[,
  `:=`(
    oo_event = fifelse(oo_t2e > 1825, 0, oo_event),
    oo_t2e = fifelse(oo_t2e > 1825, 1825, oo_t2e),
    obese_event = fifelse(obese_t2e > 1825, 0, obese_event),
    obese_t2e = fifelse(obese_t2e > 1825, 1825, obese_t2e)
  )
]

# Drop columns not needed in the analysis.
dat[,
  `:=`(
    index_date = NULL,
    state = NULL,
    visit_reason = NULL,
    ageDays = NULL,
    lenpre = NULL,
    last_visit_date = NULL,
    neuter_date_after2014 = NULL,
    dx_hypothy_date_earliest = NULL,
    dx_hyprthy_date_earliest = NULL,
    oo_date_earliest = NULL,
    obese_date_earliest = NULL,
    oo_event_date = NULL,
    obese_event_date = NULL,
    ageYearsT = NULL
  )
]

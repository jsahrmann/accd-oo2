## 03-oo-coxph.R -----------------------------------------------------
##
## Fit Cox proportional hazards models for the overweight/obese
## outcome.


### Setup ------------------------------------------------------------

library(data.table)
library(Hmisc)
library(magrittr)
library(rms)
library(survival)


### Input data -------------------------------------------------------

source("./code/00-const-fn.R")
source("./code/02-read-cohort.R")


### Data management --------------------------------------------------

##

make_table_count_all <- function(.breed) {
  dat_all_breed <- dat[
    breed == .breed,
    .N,
    by = .(breed, mixed_breed, sex, ageYearsR, sn)
  ]
  dat_all_breed <- dat_all_breed[
    order(mixed_breed, sex, ageYearsR, sn)]
  data.table::setnames(
    dat_all_breed,
    old = names(dat_all_breed),
    new = c("Breed", "MixedBreed", "Sex", "Age", "SN", "Count")
  )
}

make_table_count_oo <- function(.breed) {
  dat_oo_breed <- dat[
    breed == .breed & oo_event == 1,
    .N,
    by = .(breed, mixed_breed, sex, ageYearsR, sn)
  ]
  dat_oo_breed <- dat_oo_breed[
    order(mixed_breed, sex, ageYearsR, sn)]
  data.table::setnames(
    dat_oo_breed,
    old = names(dat_oo_breed),
    new = c("Breed", "MixedBreed", "Sex", "Age", "SN", "Count")
  )
}

make_table_count_ob <- function(.breed) {
  dat_ob_breed <- dat[
    breed == .breed & obese_event == 1,
    .N,
    by = .(breed, mixed_breed, sex, ageYearsR, sn)
  ]
  dat_ob_breed <- dat_ob_breed[
    order(mixed_breed, sex, ageYearsR, sn)]
  data.table::setnames(
    dat_ob_breed,
    old = names(dat_ob_breed),
    new = c("Breed", "MixedBreed", "Sex", "Age", "SN", "Count")
  )
}

toySmallTop10 <- sort(
  table(dat[size == "Toy and Small", breed]), decreasing = TRUE)[1:10]
largeTop5 <- sort(
  table(dat[size == "Large", breed]), decreasing = TRUE)[1:5]

counts_table_all_list <- lapply(
  names(c(toySmallTop10, largeTop5)), make_table_count_all)
names(counts_table_all_list) <- names(c(toySmallTop10, largeTop5))
counts_table_oo_list <- lapply(
  names(c(toySmallTop10, largeTop5)), make_table_count_oo)
names(counts_table_oo_list) <- names(c(toySmallTop10, largeTop5))
counts_table_ob_list <- lapply(
  names(c(toySmallTop10, largeTop5)), make_table_count_ob)
names(counts_table_ob_list) <- names(c(toySmallTop10, largeTop5))

style <- openxlsx::createStyle(halign = "center")

openxlsx::write.xlsx(
  counts_table_all_list,
  file = "./output/counts_all_MixedSexAgeSN.xlsx",
  colWidths = "auto"
)
wb_all <- openxlsx::loadWorkbook(
  "./output/counts_all_MixedSexAgeSN.xlsx")
for (i in seq_along(counts_table_all_list)) {
  nrowi <- nrow(counts_table_all_list[[i]]) + 1
  ncoli <- ncol(counts_table_all_list[[i]]) + 1
  openxlsx::addStyle(
    wb_all, sheet = i, rows = 1:nrowi, cols = 1:ncoli, style = style,
    gridExpand = TRUE
  )
}
openxlsx::saveWorkbook(
  wb_all, "./output/counts_all_MixedSexAgeSN.xlsx", overwrite = TRUE)

openxlsx::write.xlsx(
  counts_table_oo_list,
  file = "./output/counts_OO_MixedSexAgeSN.xlsx",
  colWidths = "auto"
)
wb_oo <- openxlsx::loadWorkbook(
  "./output/counts_OO_MixedSexAgeSN.xlsx")
for (i in seq_along(counts_table_oo_list)) {
  nrowi <- nrow(counts_table_oo_list[[i]]) + 1
  ncoli <- ncol(counts_table_oo_list[[i]]) + 1
  openxlsx::addStyle(
    wb_oo, sheet = i, rows = 1:nrowi, cols = 1:ncoli, style = style,
    gridExpand = TRUE
  )
}
openxlsx::saveWorkbook(
  wb_oo, "./output/counts_OO_MixedSexAgeSN.xlsx", overwrite = TRUE)

openxlsx::write.xlsx(
  counts_table_ob_list,
  file = "./output/counts_Obese_MixedSexAgeSN.xlsx",
  colWidths = "auto"
)
wb_ob <- openxlsx::loadWorkbook(
  "./output/counts_Obese_MixedSexAgeSN.xlsx")
for (i in seq_along(counts_table_ob_list)) {
  nrowi <- nrow(counts_table_ob_list[[i]]) + 1
  ncoli <- ncol(counts_table_ob_list[[i]]) + 1
  openxlsx::addStyle(
    wb_ob, sheet = i, rows = 1:nrowi, cols = 1:ncoli, style = style,
    gridExpand = TRUE
  )
}
openxlsx::saveWorkbook(
  wb_ob, "./output/counts_Obese_MixedSexAgeSN.xlsx", overwrite = TRUE)

#

xtabs(
  ~ ageYearsR + sn + sex,
  data = dat, subset = breed == "Pug" & oo_event == 1
)
xtabs(
  ~ ageYearsR + sn + sex,
  data = dat, subset = breed == "Pug" & obese_event == 1
)

xtabs(
  ~ ageYearsR + sn + sex + mixed_breed,
  data = dat, subset = breed == "Pug" & oo_event == 1
)
x <- xtabs(
  ~ ageYearsR + sn + sex + mixed_breed,
  data = dat, subset = breed == "Pug" & obese_event == 1
)

##

toySmallTop10 <- sort(
  table(dat[size == "Toy and Small", breed]), decreasing = TRUE)[1:10]
largeTop5 <- sort(
  table(dat[size == "Large", breed]), decreasing = TRUE)[1:5]

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
  breed_dat <- dat[breed == .breed]
  sans_breed_dat <- dat[breed != .breed]
  list(
    dat = .data,
    breed_dat = breed_dat,
    sans_breed_dat = sans_breed_dat)
}

oo2_fit_models <- function(ds_list) {
  # Define data distributions for each data set.  We need to super
  # assign these because =rms= isn't able to find them when we run
  # ~options(datadist = "dd")~ unless they're in the global
  # environment.
  dat_datadist <<- rms::datadist(ds_list$dat)
  breed_datadist <<- rms::datadist(ds_list$breed_dat)
  sans_breed_datadist <<- rms::datadist(ds_list$sans_breed_dat)
  # Fit the O/O model for each data set.
  dat_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$dat, x = TRUE, y = TRUE, surv = TRUE
  )
  breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  sans_breed_oo_mod <- rms::cph(
    survival::Surv(oo_t2e, oo_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$sans_breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
  # Fit the obese model for each data set.
  dat_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$dat, x = TRUE, y = TRUE, surv = TRUE
  )
  breed_ob_mod <- rms::cph(
    survival::Surv(obese_t2e, obese_event) ~
      sn + rcs(ageYearsX, 5) + sex + mixed_breed + wellness_plan +
      rcs(weight, 3) + rcs(visitsPerYear, 3) +
      sn : rcs(ageYearsX, 5) + sn : sex + sn : rcs(weight, 3) +
      sex : rcs(ageYearsX, 5) + sex : rcs(weight, 3),
    data = ds_list$breed_dat, x = TRUE, y = TRUE, surv = TRUE
  )
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
    dat = list(
      oo_mod = dat_oo_mod,
      ob_mod = dat_ob_mod
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

oo2_sn_evaluate_models <- function(model_list, size_category, breed) {
  # Create the data set of reference points at which to evaluate the
  # models.
  ref_pts <- data.table::as.data.table(
    define_sn_reference_points(
      sizes = size_category,
      ages = c(0.25, seq(0.5, 6, by = 0.5))
    )
  )
  # Evaluate the models.
  options(datadist = "dat_datadist")
  dat_oo_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$dat$oo_mod)
  dat_ob_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$dat$ob_mod)
  options(datadist = "breed_datadist")
  breed_oo_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$breed$oo_mod)
  breed_ob_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$breed$ob_mod)
  options(datadist = "sans_breed_datadist")
  sans_breed_oo_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$sans_breed$oo_mod)
  sans_breed_ob_est <- evaluate_sn_reference_points(
    ref_pts, model = model_list$sans_breed$ob_mod)
  # Combine the reference points and model estimates into individual
  # data sets.
  dat_oo <- cbind(ref_pts, dat_oo_est)[, analysis := size_category]
  breed_oo <- cbind(ref_pts, breed_oo_est)[, analysis := breed]
  sans_breed_oo <- cbind(ref_pts, sans_breed_oo_est)[,
    analysis := paste0(size_category, ", w/o ", breed)
  ]
  dat_ob <- cbind(ref_pts, dat_ob_est)[, analysis := size_category]
  breed_ob <- cbind(ref_pts, breed_ob_est)[, analysis := breed]
  sans_breed_ob <- cbind(ref_pts, sans_breed_ob_est)[,
    analysis := paste0(size_category, ", w/o ", breed)
  ]
  # Stack the outcome-specific data sets for export.
  list(
    oo_effect_data = rbind(dat_oo, breed_oo, sans_breed_oo),
    ob_effect_data = rbind(dat_ob, breed_ob, sans_breed_ob)
  )
}

oo2_sn_plot_effects <- function(
  effect_data_list, size_category, breed, pt_offset = 0.05,
  line_size = 0.4, rel_point_size = 1.5
) {
  # Subset to just estimates at median weight, shift age so that
  # points in the effect plot don't overlap, and convert ~analysis~ to
  # a factor to control ordering in the legend.
  oo_effect_data_for_plot <- effect_data_list$oo_effect_data[
    wt_pctl == 50
  ][,
    `:=`(
      age = data.table::fcase(
        analysis == size_category,
        age,
        analysis == paste0(size_category, ", w/o ", breed),
        age - pt_offset,
        analysis == breed,
        age + pt_offset
      ),
      analysis = factor(
        analysis,
        levels = c(
          size_category, breed, paste0(size_category, ", w/o ", breed)
        )
      )
    )
  ]
  ob_effect_data_for_plot <- effect_data_list$ob_effect_data[
    wt_pctl == 50
  ][,
    `:=`(
      age = data.table::fcase(
        analysis == size_category,
        age,
        analysis == paste0(size_category, ", w/o ", breed),
        age - pt_offset,
        analysis == breed,
        age + pt_offset
      ),
      analysis = factor(
        analysis,
        levels = c(
          size_category, breed, paste0(size_category, ", w/o ", breed)
        )
      )
    )
  ]
  # Create the effect plots.
  oo_plot <- ggplot(
    data = oo_effect_data_for_plot,
    mapping = aes(
      x = age, y = hr, ymin = lo, ymax = hi, colour = analysis)
  ) +
    geom_pointrange(size = line_size, fatten = rel_point_size) +
    geom_line() +
    scale_y_log10() +
    scale_colour_discrete("") +
    ggtitle(
      paste("Effect of S/N on O/O Outcome,", breed, "Analysis")
    ) +
    xlab("\nAge at Index (Years)") +
    ylab("Hazard Ratio for Gonadectomy\n") +
    facet_wrap(vars(sex)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
  ob_plot <- ggplot(
    data = ob_effect_data_for_plot,
    mapping = aes(
      x = age, y = hr, ymin = lo, ymax = hi, colour = analysis)
  ) +
    geom_pointrange(size = line_size, fatten = rel_point_size) +
    geom_line() +
    scale_y_log10() +
    scale_colour_discrete("") +
    ggtitle(
      paste("Effect of S/N on Obese-Only Outcome,", breed, "Analysis")
    ) +
    xlab("\nAge at Index (Years)") +
    ylab("Hazard Ratio for Gonadectomy\n") +
    facet_wrap(vars(sex)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
  list(
    oo_plot = oo_plot,
    ob_plot = ob_plot
  )
}

oo2_age_at_sn_evaluate_models <- function(
  model_list, size_category, breed
) {
  # Create the data set of reference points at which to evaluate the
  # models.
  ref_pts <- data.table::as.data.table(
    define_age_reference_points(
      sizes = size_category,
      ages = c(0.25, seq(0.5, 6, by = 0.5))
    )
  )
  # Evaluate the models.
  options(datadist = "dat_datadist")
  dat_oo_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$dat$oo_mod)
  dat_ob_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$dat$ob_mod)
  options(datadist = "breed_datadist")
  breed_oo_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$breed$oo_mod)
  breed_ob_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$breed$ob_mod)
  options(datadist = "sans_breed_datadist")
  sans_breed_oo_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$sans_breed$oo_mod)
  sans_breed_ob_est <- evaluate_age_among_sn_reference_points(
    ref_pts, model = model_list$sans_breed$ob_mod)
  # Combine the reference points and model estimates into individual
  # data sets.
  dat_oo <- cbind(ref_pts, dat_oo_est)[, analysis := size_category]
  breed_oo <- cbind(ref_pts, breed_oo_est)[, analysis := breed]
  sans_breed_oo <- cbind(ref_pts, sans_breed_oo_est)[,
    analysis := paste0(size_category, ", w/o ", breed)
  ]
  dat_ob <- cbind(ref_pts, dat_ob_est)[, analysis := size_category]
  breed_ob <- cbind(ref_pts, breed_ob_est)[, analysis := breed]
  sans_breed_ob <- cbind(ref_pts, sans_breed_ob_est)[,
    analysis := paste0(size_category, ", w/o ", breed)
  ]
  # Stack the outcome-specific data sets for export.
  list(
    oo_effect_data = rbind(dat_oo, breed_oo, sans_breed_oo),
    ob_effect_data = rbind(dat_ob, breed_ob, sans_breed_ob)
  )
}

oo2_age_at_sn_plot_effects <- function(
  effect_data_list, size_category, breed, pt_offset = 0.05,
  line_size = 0.4, rel_point_size = 1.5
) {
  # Subset to just estimates with a reference age of one year, shift
  # comparator age so that points in the effect plot don't overlap,
  # and convert ~analysis~ to a factor to control ordering in the
  # legend.
  oo_effect_data_for_plot <- effect_data_list$oo_effect_data[
    reference_age == 1.00
  ][,
    `:=`(
      comparator_age = data.table::fcase(
        analysis == size_category,
        comparator_age,
        analysis == paste0(size_category, ", w/o ", breed),
        comparator_age - pt_offset,
        analysis == breed,
        comparator_age + pt_offset
      ),
      analysis = factor(
        analysis,
        levels = c(
          size_category, breed, paste0(size_category, ", w/o ", breed)
        )
      )
    )
  ]
  ob_effect_data_for_plot <- effect_data_list$ob_effect_data[
    reference_age == 1.00
  ][,
    `:=`(
      comparator_age = data.table::fcase(
        analysis == size_category,
        comparator_age,
        analysis == paste0(size_category, ", w/o ", breed),
        comparator_age - pt_offset,
        analysis == breed,
        comparator_age + pt_offset
      ),
      analysis = factor(
        analysis,
        levels = c(
          size_category, breed, paste0(size_category, ", w/o ", breed)
        )
      )
    )
  ]
  # Create the effect plots.
  oo_plot <- ggplot(
    data = oo_effect_data_for_plot,
    mapping = aes(
      x = comparator_age, y = hr, ymin = lo, ymax = hi,
      colour = analysis
    )
  ) +
    geom_pointrange(size = line_size, fatten = rel_point_size) +
    geom_line() +
    scale_y_log10() +
    scale_colour_discrete("") +
    ggtitle(
      paste("Effect of Age at S/N on O/O Outcome,", breed, "Analysis")
    ) +
    xlab("\nAge at Gonadectomy (Years)") +
    ylab("Hazard Ratio for Age\n") +
    facet_wrap(vars(sex)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
  ob_plot <- ggplot(
    data = ob_effect_data_for_plot,
    mapping = aes(
      x = comparator_age, y = hr, ymin = lo, ymax = hi,
      colour = analysis
    )
  ) +
    geom_pointrange(size = line_size, fatten = rel_point_size) +
    geom_line() +
    scale_y_log10() +
    scale_colour_discrete("") +
    ggtitle(
      paste(
        "Effect of Age at S/N on Obese-Only Outcome,", breed,
        "Analysis"
      )
    ) +
    xlab("\nAge at Gonadectomry (Years)") +
    ylab("Hazard Ratio for Age\n") +
    facet_wrap(vars(sex)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
  list(
    oo_plot = oo_plot,
    ob_plot = ob_plot
  )
}


#### Toy and Small, Pug ----------------------------------------------

pug_ds <- oo2_make_data_sets(dat, "Toy and Small", "Pug")
pug_models <- oo2_fit_models(pug_ds)
pug_sn_effect_data <- oo2_sn_evaluate_models(
  pug_models, size_category = "Toy and Small", breed = "Pug")
pug_sn_effect_plots <- oo2_sn_plot_effects(
  pug_sn_effect_data,
  size_category = "Toy and Small", breed = "Pug")
pug_age_at_sn_effect_data <- oo2_age_at_sn_evaluate_models(
  pug_models, size_category = "Toy and Small", breed = "Pug")
pug_age_at_sn_effect_plots <- oo2_age_at_sn_plot_effects(
  pug_age_at_sn_effect_data,
  size_category = "Toy and Small", breed = "Pug"
)

cairo_pdf(
  "./output/fig/20230701_ToySmall_Pug_1_SN_OO_Figure.pdf",
  width = 8, height = 6
)
pug_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/20230701_ToySmall_Pug_2_SN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
pug_sn_effect_plots$ob_plot
dev.off()

cairo_pdf(
  "./output/fig/20230701_ToySmall_Pug_3_AgeAtSN_OO_Figure.pdf",
  width = 8, height = 6
)
pug_age_at_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/20230701_ToySmall_Pug_4_AgeAtSN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
pug_age_at_sn_effect_plots$ob_plot
dev.off()


#### Toy and Small, Yorkshire Terrier --------------------------------

yorkt_ds <- oo2_make_data_sets(
  dat, "Toy and Small", "Yorkshire Terrier")
yorkt_models <- oo2_fit_models(yorkt_ds)
yorkt_sn_effect_data <- oo2_sn_evaluate_models(
  yorkt_models, "Toy and Small", "Yorkshire Terrier")
yorkt_sn_effect_plots <- oo2_sn_plot_effects(
  yorkt_sn_effect_data, "Toy and Small", "Yorkshire Terrier")
yorkt_age_at_sn_effect_data <- oo2_age_at_sn_evaluate_models(
  yorkt_models, "Toy and Small", "Yorkshire Terrier")
yorkt_age_at_sn_effect_plots <- oo2_age_at_sn_plot_effects(
  yorkt_age_at_sn_effect_data, "Toy and Small", "Yorkshire Terrier")


cairo_pdf(
  "./output/fig/ToySmall_YorkshireTerrier_1_SN_OO_Figure.pdf",
  width = 8, height = 6
)
yorkt_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_YorkshireTerrier_2_SN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
yorkt_sn_effect_plots$ob_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_YorkshireTerrier_3_AgeAtSN_OO_Figure.pdf",
  width = 8, height = 6
)
yorkt_age_at_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_YorkshireTerrier_4_AgeAtSN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
yorkt_age_at_sn_effect_plots$ob_plot
dev.off()


#### Toy and Small, Chihuahua ----------------------------------------

chih_ds <- oo2_make_data_sets(
  dat, "Toy and Small", "Chihuahua")
chih_models <- oo2_fit_models(chih_ds)
chih_sn_effect_data <- oo2_sn_evaluate_models(
  chih_models, "Toy and Small", "Chihuahua")
chih_sn_effect_plots <- oo2_sn_plot_effects(
  chih_sn_effect_data, "Toy and Small", "Chihuahua")
chih_age_at_sn_effect_data <- oo2_age_at_sn_evaluate_models(
  chih_models, "Toy and Small", "Chihuahua")
chih_age_at_sn_effect_plots <- oo2_age_at_sn_plot_effects(
  chih_age_at_sn_effect_data, "Toy and Small", "Chihuahua")


cairo_pdf(
  "./output/fig/ToySmall_Chihuahua_1_SN_OO_Figure.pdf",
  width = 8, height = 6
)
chih_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_Chihuahua_2_SN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
chih_sn_effect_plots$ob_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_Chihuahua_3_AgeAtSN_OO_Figure.pdf",
  width = 8, height = 6
)
chih_age_at_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/ToySmall_Chihuahua_4_AgeAtSN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
chih_age_at_sn_effect_plots$ob_plot
dev.off()


#### Large, American Bulldog -----------------------------------------

ambd_ds <- oo2_make_data_sets(
  dat, "Large", "American Bulldog")
ambd_models <- oo2_fit_models(ambd_ds)
ambd_sn_effect_data <- oo2_sn_evaluate_models(
  ambd_models, "Large", "American Bulldog")
ambd_sn_effect_plots <- oo2_sn_plot_effects(
  ambd_sn_effect_data, "Large", "American Bulldog")
ambd_age_at_sn_effect_data <- oo2_age_at_sn_evaluate_models(
  ambd_models, "Large", "American Bulldog")
ambd_age_at_sn_effect_plots <- oo2_age_at_sn_plot_effects(
  ambd_age_at_sn_effect_data, "Large", "American Bulldog")


cairo_pdf(
  "./output/fig/Large_AmericanBulldog_1_SN_OO_Figure.pdf",
  width = 8, height = 6
)
ambd_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_AmericanBulldog_2_SN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
ambd_sn_effect_plots$ob_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_AmericanBulldog_3_AgeAtSN_OO_Figure.pdf",
  width = 8, height = 6
)
ambd_age_at_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_AmericanBulldog_4_AgeAtSN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
ambd_age_at_sn_effect_plots$ob_plot
dev.off()


#### Large, Labrador Retriever -----------------------------------------

labr_ds <- oo2_make_data_sets(
  dat, "Large", "Labrador Retriever")
labr_models <- oo2_fit_models(labr_ds)
labr_sn_effect_data <- oo2_sn_evaluate_models(
  labr_models, "Large", "Labrador Retriever")
labr_sn_effect_plots <- oo2_sn_plot_effects(
  labr_sn_effect_data, "Large", "Labrador Retriever")
labr_age_at_sn_effect_data <- oo2_age_at_sn_evaluate_models(
  labr_models, "Large", "Labrador Retriever")
labr_age_at_sn_effect_plots <- oo2_age_at_sn_plot_effects(
  labr_age_at_sn_effect_data, "Large", "Labrador Retriever")


cairo_pdf(
  "./output/fig/Large_LabradorRetriever_1_SN_OO_Figure.pdf",
  width = 8, height = 6
)
labr_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_LabradorRetriever_2_SN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
labr_sn_effect_plots$ob_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_LabradorRetriever_3_AgeAtSN_OO_Figure.pdf",
  width = 8, height = 6
)
labr_age_at_sn_effect_plots$oo_plot
dev.off()

cairo_pdf(
  "./output/fig/Large_LabradorRetriever_4_AgeAtSN_ObeseOnly_Figure.pdf",
  width = 8, height = 6
)
labr_age_at_sn_effect_plots$ob_plot
dev.off()

pdftools::pdf_combine(
  input = dir("./output/fig", "LabradorRetriever", full.names = TRUE),
  output = "./output/fig/Large_LabradorRetriever_Figures.pdf"
)
pdftools::pdf_combine(
  input = dir("./output/fig", "AmericanBulldog", full.names = TRUE),
  output = "./output/fig/Large_AmericanBulldog_Figures.pdf"
)
pdftools::pdf_combine(
  input = dir("./output/fig", "YorkshireTerrier", full.names = TRUE),
  output = "./output/fig/ToySmall_YorkshireTerrier_Figures.pdf"
)
pdftools::pdf_combine(
  input = dir("./output/fig", "Pug", full.names = TRUE),
  output = "./output/fig/ToySmall_Pug_Figures.pdf"
)

# Manipulates hanam data

library(tidyverse)

# Directories to be used later
data_dir <- here::here("data")
data_raw_dir <- here::here("data-raw")

# Functions ===================================================================

read_hanam_raw <- function(name) {
  read_csv(
    file.path(data_raw_dir, paste0(name, ".csv")),
    col_types = cols_only(
      hhold = col_character(),
      ind = col_character(),
      season = col_integer(),
      virus = col_character(),
      status = col_character(),
      preHI = col_integer(),
      dob = col_date(),
      measure_date = col_date()
    )
  ) %>%
    rename(prehi = preHI)
}

save_hanam <- function(dat, name) {
  write_csv(dat, file.path(data_dir, paste0("hanam-", name, ".csv")))
}

# Script ======================================================================

hanam <- read_hanam_raw("hanam")

hanam_extra <- hanam %>%
  mutate(
    status_bin = if_else(status == "Not infected", 0, 1),
    status_bin_lbl = factor(
      status_bin,
      levels = c(0, 1), labels = c("Not infected", "Infected")
    ),
    loghi = log(prehi),
    loghilb = if_else(prehi == 5L | is.na(prehi), -1e6, loghi),
    loghiub = if_else(prehi == 1280L | is.na(prehi), 1e6, loghi + log(2)),
    loghimid = case_when(
      is.na(prehi) ~ NA_real_,
      prehi == 5L ~ log(5),
      prehi == 1280L ~ log(1280),
      TRUE ~ loghi + log(2) / 2
    ),
    age_years = (measure_date - dob) / 365.25
  )

hanam_hi_info <- filter(hanam_extra, !is.na(prehi) | !is.na(status))

hanam_hi_general <- mutate(hanam_hi_info, population = "General")
save_hanam(hanam_hi_general, "hi-gen")

hanam_hi_exposed <- hanam_hi_info %>%
  group_by(virus, hhold, season) %>%
  filter(sum(status != "Not infected", na.rm = TRUE) > 0) %>%
  mutate(population = "Exposed")
save_hanam(hanam_hi_exposed, "hi-exp")

hanam_hi_summ <- bind_rows(hanam_hi_general, hanam_hi_exposed) %>%
  filter(!is.na(prehi), !is.na(status)) %>%
  group_by(population, virus, prehi, loghimid) %>%
  summarise(ntot = n(), inf_prop = sum(status != "Not infected") / ntot) %>%
  ungroup()
save_hanam(hanam_hi_summ, "hi-summ")

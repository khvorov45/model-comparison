# Manipulates hanam data
# Arseniy Khvorov
# Created 2019/12/11
# Last edit 2020/01/28

library(tidyverse)

# Directories to be used later
data_dir <- "data"
data_raw_dir <- "data-raw"

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
      preHI = col_integer()
    )
  )
}

save_hanam <- function(dat, name) {
  write_csv(dat, file.path(data_dir, paste0("hanam-", name, ".csv")))
}

# Script ======================================================================

hanam <- read_hanam_raw("hanam")

hanam_extra <- hanam %>%
  mutate(
    status_bin = if_else(status == "Not infected", 0, 1),
    logHI = log(preHI),
    logHIlb = if_else(preHI == 5L | is.na(preHI), -1e6, logHI),
    logHIub = if_else(preHI == 1280L | is.na(preHI), 1e6, logHI + log(2)),
    logHImid = case_when(
      is.na(preHI) ~ NA_real_,
      preHI == 5L ~ log(5),
      preHI == 1280L ~ log(1280),
      TRUE ~ logHI + log(2) / 2
    )
  )

hanam_hi_info <- filter(hanam_extra, !is.na(preHI) | !is.na(status))

hanam_hi_general <- mutate(hanam_hi_info, population = "General")
save_hanam(hanam_hi_general, "hi-gen")

hanam_hi_exposed <- hanam_hi_info %>%
  group_by(virus, hhold, season) %>%
  filter(sum(status != "Not infected", na.rm = TRUE) > 0) %>%
  mutate(population = "Exposed")
save_hanam(hanam_hi_exposed, "hi-exp")

hanam_hi_summ <- bind_rows(hanam_hi_general, hanam_hi_exposed) %>%
  filter(!is.na(preHI), !is.na(status)) %>%
  group_by(population, virus, preHI, logHImid) %>%
  summarise(ntot = n(), inf_prop = sum(status != "Not infected") / ntot) %>%
  ungroup()
save_hanam(hanam_hi_summ, "hi-summ")

# Manipulates hanam data
# Arseniy Khvorov
# Created 2019/12/11
# Last edit 2019/12/11

library(tidyverse)

# Directories to be used later
data_dir <- "data"

# Functions ===================================================================

save_hanam <- function(dat, name) {
  write_csv(dat, file.path(data_dir, paste0("hanam-", name, ".csv")))
}

# Script ======================================================================

hanam <- read_csv(file.path(data_dir, "hanam.csv"))

hanam_extra <- hanam %>%
  mutate(
    status_bin = if_else(status == "Not infected", 0, 1),
    logHIlb = if_else(near(preHI, 5) | is.na(preHI), -1e6, log(preHI)),
    logHIub = if_else(near(preHI, 1280) | is.na(preHI), 1e6, log(2 * preHI))
  )

hanam_hi_info <- filter(hanam_extra, !is.na(preHI) | !is.na(status))

hanam_hi_general <- mutate(hanam_hi_info, population = "general")
save_hanam(hanam_hi_general, "HI-gen")

hanam_hi_exposed <- hanam_hi_info %>%
  group_by(virus, hhold, season) %>%
  filter(sum(status != "Not infected", na.rm = TRUE) > 0) %>%
  mutate(population = "exposed")
save_hanam(hanam_hi_exposed, "HI-exp")

hanam_hi_summ <- bind_rows(hanam_hi_general, hanam_hi_exposed) %>%
  filter(!is.na(preHI), !is.na(status)) %>%
  group_by(population, virus, preHI) %>%
  summarise(ntot = n(), inf_prop = sum(status != "Not infected") / ntot)
save_hanam(hanam_hi_summ, "HI-summ")

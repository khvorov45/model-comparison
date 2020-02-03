# Relevant graphs of hanam data
# Arseniy Khvorov
# Created 2020-02-03
# Last edit 2020-02-03

library(tidyverse)

# Folders to be used later
data_dir <- "data"
data_table_dir <- "data-table"

# Functions ===================================================================

read_hanam <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")), 
    col_types = cols(
      hhold = col_character(),
      ind = col_character(),
      season = col_integer(),
      virus = col_character(),
      status = col_character(),
      prehi = col_integer(),
      status_bin = col_integer(),
      loghi = col_double(),
      loghilb = col_double(),
      loghimid = col_double(),
      population = col_factor(levels = c("General", "Exposed")),
      age_years = col_double()
    )
  )
}

save_table <- function(dat, name) {
  write_csv(dat, file.path(data_table_dir, paste0(name, ".csv")))
}

# Script ======================================================================

han_dat <- map_dfr(c("hanam-hi-gen", "hanam-hi-exp"), read_hanam) %>%
  filter(!is.na(status), !is.na(prehi), virus != "H1N1seas")

han_tbl1 <- han_dat %>%
  group_by(population, virus, prehi) %>%
  summarise(
    prot0 = sum(status_bin == 0L),
    prop0 = prot0 / n(),
    mean.age0 = mean(age_years[status_bin == 0L], na.rm = TRUE),
    prot1 = sum(status_bin),
    prop1 = prot1 / n(),
    mean.age1 = mean(age_years[status_bin == 1L], na.rm = TRUE)
  )

save_table(han_tbl1, "hanam-hi-tbl1")

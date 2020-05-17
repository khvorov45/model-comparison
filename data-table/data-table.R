# Relevant tables of data

library(tidyverse)

# Folders to be used later
data_dir <- here::here("data")
data_table_dir <- here::here("data-table")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

save_table <- function(dat, name) {
  write_csv(dat, file.path(data_table_dir, paste0(name, ".csv")))
}

# Script ======================================================================

han_dat <- map_dfr(c("hanam-hi-gen", "hanam-hi-exp"), read_data) %>%
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

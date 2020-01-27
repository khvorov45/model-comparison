# Kiddyvax data manipulation
# Arseniy Khvorov
# Created 2020-01-17
# Last edit 2020-01-28

library(tidyverse)

# Directories to be used later
data_raw_dir <- "data-raw"
data_dir <- "data"

# Functions ===================================================================

read_swab <- function(filepath) {
  filepath %>%
    read_csv(
      col_types = cols_only(
        hhID = col_integer(), # one observation per household
        date = col_date("%Y%m%d"),
        FluA = col_character(),
        `Swine H1` = col_character(),
        H1 = col_character(),
        H3 = col_character(),
        FluB = col_character(),
        `FluB subtype` = col_character()
      )
    ) %>%
    select(
      id = hhID, swab_date = date, a = FluA,
      h1pdm = `Swine H1`, h1seas = H1, h3 = H3, b = FluB,
      b_type = `FluB subtype`
    )
}

fix_subtypes_swab <- function(swab) {
  swab %>%
    mutate(
      bvic = case_when(
        b == "N" ~ "N",
        b_type == "Victoria" ~ "P",
        b_type == "Yamagata" ~ "N",
        is.na(b) ~ NA_character_
      ),
      byam = case_when(
        b == "N" ~ "N",
        b_type == "Victoria" ~ "N",
        b_type == "Yamagata" ~ "Y",
        is.na(b) ~ NA_character_
      ),
      h1pdm = if_else(a == "N", "N", h1pdm),
      h1seas = if_else(a == "N", "N", h1seas),
      h3 = if_else(a == "N", "N", h3),
    )
}

lengthen_swab <- function(swab) {
  swab %>%
    select(-a, -b_type, -b) %>%
    pivot_longer(
      cols = c(h1pdm, h1seas, h3, bvic, byam),
      names_to = "virus", values_to = "swab_result"
    ) %>%
    mutate(swab_result = if_else(swab_result == "N", 0L, 1L))
}

# Tests that the reshaped data contains the same information
test_swab <- function(filepath) {
  swab_og <- read_swab(filepath)
  swab <- swab_og %>% fix_a_subtypes_swab() %>% lengthen_swab()
  
  h1pdm_og <- swab_og %>% filter(h1pdm == "P") %>% pull(id)
  h1pdm <- swab %>% filter(swab_result == 1 & virus == "h1pdm") %>% pull(id)
  testthat::expect_equal(h1pdm_og, h1pdm)
  
  h1seas_og <- swab_og %>% filter(h1seas == "P") %>% pull(id)
  h1seas <- swab %>% filter(swab_result == 1 & virus == "h1seas") %>% pull(id)
  testthat::expect_equal(h1seas_og, h1seas)
  
  h3_og <- swab_og %>% filter(h3 == "P") %>% pull(id)
  h3 <- swab %>% filter(swab_result == 1 & virus == "h3") %>% pull(id)
  testthat::expect_equal(h3_og, h3)
  
  bvic_og <- swab_og %>% filter(b_type == "Victoria") %>% pull(id)
  bvic <- swab %>% filter(swab_result == 1 & virus == "bvic") %>% pull(id)
  testthat::expect_equal(bvic_og, bvic)
  
  byam_og <- swab_og %>% filter(b_type == "Yamagata") %>% pull(id)
  byam <- swab %>% filter(swab_result == 1 & virus == "byam") %>% pull(id)
  testthat::expect_equal(byam_og, byam)
}

# Say they were infected if any of the swabs are positive
# i.e. status is infection at any time during the study
summarise_swab <- function(swab) {
  find_inf_date <- function(swab_result, swab_date) {
    if (all(swab_result == 0)) return(as.Date(NA))
    first(swab_date[swab_result == 1])
  }
  swab %>%
    filter(swab_date <= end_date, !is.na(swab_result)) %>%
    group_by(id, virus) %>%
    summarise(
      status = as.integer(any(swab_result == 1)),
      infection_date = find_inf_date(swab_result, swab_date)
    ) %>%
    ungroup()
}

read_serology <- function(filepath) {
  filepath %>%
    read_csv(
      col_types = cols_only(
        hhID = col_integer(),
        start.date = col_date("%d/%m/%Y"),
        end.date = col_date("%d/%m/%Y"),
        postvax.pH1 = col_integer(),
        postvax.sH1 = col_integer(),
        postvax.sH3 = col_integer(),
        postvax.B.Brisbane = col_integer(),
        postvax.B.Floride = col_integer()
      )
    ) %>%
    rename(id = hhID, start_date = start.date, end_date = end.date)
}

lengthen_serology <- function(serology) {
  serology %>%
    pivot_longer(
      starts_with("postvax."), 
      names_to = "virus", values_to = "hi"
    )
}

fix_subtypes_serology <- function(serology) {
  serology %>%
    mutate(
      virus = str_replace(virus, "postvax.", "") %>% 
        recode(
          "pH1" = "h1pdm", "sH1" = "h1seas", "sH3" = "h3",
          "B.Brisbane" = "bvic", "B.Floride" = "byam"
        )
    )
}

add_loghis <- function(dat) {
  dat %>%
    mutate(
      loghi = log(hi),
      loghimid = if_else(hi == 5, log(5), loghi + log(2) / 2),
      loghilb = if_else(hi == 5 | is.na(hi), -1e6, loghi),
      loghiub = if_else(is.na(hi), 1e6, loghi + log(2))
    )
}

save_data <- function(dat, name) {
  write_csv(dat, file.path(data_dir, paste0(name, ".csv")))
}

# Script ======================================================================

swab <- read_swab(file.path(data_raw_dir, "kiddyvax-swab.csv")) %>% 
  fix_subtypes_swab() %>%
  lengthen_swab()

serology <- read_serology(file.path(data_raw_dir, "kiddyvax-serology.csv")) %>% 
  lengthen_serology() %>%
  fix_subtypes_serology() %>%
  add_loghis()

startend_dates <- select(serology, id, start_date, end_date) %>% unique()

swab_full <- left_join(swab, startend_dates, by = "id")
save_data(swab_full, "kiddyvax-swab")

swab_summ <- summarise_swab(swab_full)

full_data <- full_join(swab_summ, serology, by = c("id", "virus"))
save_data(full_data, "kiddyvax-main")

full_data_summ <- full_data %>%
  group_by(virus, hi, loghimid) %>%
  filter(!is.na(status), !is.na(hi)) %>%
  summarise(ntot = n(), inf_prop = sum(status) / ntot) %>%
  ungroup()
save_data(full_data_summ, "kiddyvax-main-summ")

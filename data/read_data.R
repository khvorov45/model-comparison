read_data <- function(name) {
  all_col_types <- list(
    hanam = cols(
      hhold = col_character(),
      ind = col_character(),
      season = col_integer(),
      virus = col_character(),
      status = col_character(),
      prehi = col_integer(),
      status_bin = col_integer(),
      status_bin_lbl = col_factor(c("Not infected", "Infected")),
      loghi = col_double(),
      loghilb = col_double(),
      loghimid = col_double(),
      population = col_factor(levels = c("General", "Exposed"))
    ),
    hanam_summ = cols(
      population = col_factor(levels = c("General", "Exposed")),
      virus = col_character(),
      prehi = col_integer(),
      loghimid = col_double(),
      ntot = col_integer(),
      inf_prop = col_double()
    )
  )
  if (name == "hanam-hi-exp" | name == "hanam-hi-gen") {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = all_col_types[["hanam"]]
    )
  } else if (name == "hanam-hi-summ") {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = all_col_types[["hanam_summ"]]
    )
  } else {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = cols()
    )
  }
}

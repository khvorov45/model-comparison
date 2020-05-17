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
    ),
    kiddyvaxmain = cols(
      id = col_integer(),
      virus = col_character(),
      status = col_integer(),
      hi = col_integer(),
      virus_lbl = col_factor(c("B Vic", "B Yam", "A H1pdm", "A H1seas", "A H3"))
    ),
    swab = cols(
      id = col_integer(),
      swab_date = col_date("%Y-%m-%d"),
      virus = col_character(),
      swab_result = col_integer(),
      start_date = col_date("%Y-%m-%d"),
      end_date = col_date("%Y-%m-%d"),
      virus_lbl = col_factor(c("B Vic", "B Yam", "A H1pdm", "A H1seas", "A H3"))
    ),
    kv_summ = cols(
      hi = col_integer(),
      ntot = col_integer(),
      virus_lbl = col_factor(c("B Vic", "B Yam", "A H1pdm", "A H1seas", "A H3"))
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
  } else if (name == "kiddyvaxmain") {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = all_col_types[["kiddyvaxmain"]]
    )
  } else if (name == "kiddyvaxmain-swab") {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = all_col_types[["swab"]]
    )
  } else if (name == "kiddyvaxmain-summ") {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = all_col_types[["kv_summ"]]
    )
  } else {
    read_csv(
      file.path("data", glue::glue("{name}.csv")),
      col_types = cols()
    )
  }
}

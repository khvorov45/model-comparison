# Extract the data behind Cox fits of the lower panel of Figure 1, Ng (2013)
# Arseniy Khvorov
# Created 2020-02-20
# Last edit 2020-02-20

library(tidyverse)

# Directories used
data_dir <- "data"

# Functions ===================================================================

save_data <- function(data, name) {
  write_csv(data, file.path(data_dir, glue::glue("{name}.csv")))
}

mod_soph_dat <- function(data, model_name) {
  names(data) <- tolower(names(data))
  postvax_name <- names(data)[str_detect(names(data), "postvax")]
  as_tibble(data) %>%
    group_by(hhid) %>%
    mutate(
      start = t,
      end = lead(t, default = max(t) + 1),
      model = model_name
    ) %>%
    ungroup() %>%
    rename(postvax = !!rlang::sym(postvax_name))
}

# Sophia's script =============================================================

basedir <- "http://web.hku.hk/~bcowling/influenza"
scriptdir <- file.path(basedir, "KiddivaxMain_JID_titerprot_scripts")
datadir <- file.path(basedir, "data")
dir <- file.path(datadir, "KiddivaxMainV2")

source(file.path(scriptdir, "source_1.r"))
source(file.path(scriptdir, "source_2.r"))

library(survival)

#
# Pandemic A(H1N1)  Model
#

# Read in data
kdata <- read.csv(file.path(dir, "serology.csv"))
swab <- read.csv(file.path(dir, "swab.csv"))
random <- read.csv(file.path(dir, "randomcode.csv"))
fluproxy <- read.csv(file.path(datadir, "influenza_proxy_1998to2013.csv"))

time0 <- "1/10/2009"
fluproxy <- fluproxy[614:679, ]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending, format = "%d/%m/%Y") - as.Date(time0, format = "%d/%m/%Y")) - 1
proxy <- data.frame(ti = fluproxy$ti, proxy = (1 / r) * fluproxy$A.H1N1pdm.proxy)

# create smoothed proxy for background risk
s.proxy <- data.frame(
  ti = 1:max(proxy$ti), 
  s.proxy = predict(
    smooth.spline(proxy$ti, proxy$proxy, spar = 0.6), 1:max(proxy$ti)
  )$y
)

kdata <- kdata[kdata$member == 0, ]
kdata <- merge(kdata, random, by = "hhID", all.x = T)

swab <- swab[swab$Swine.H1 == "P", ]
swab <- swab[order(swab$hhID, swab$member, swab$date), ]
swab$mark <- 0
for (i in 2:nrow(swab)) {
  if (swab$hhID[i] == swab$hhID[i - 1] & swab$member[i] == swab$member[i - 1]) swab$mark[i] <- 1
}
swab <- swab[swab$mark == 0 & swab$member == 0, ]

# merge in PCR type and date
kdata$swab.pH1 <- 1 * (kdata$hhID %in% swab$hhID[swab$Swine.H1 == "P"])
kdata <- kdata[!is.na(kdata$postvax.pH1) & kdata$postvax != "", ]
kdata2 <- merge(kdata, swab[swab$Swine.H1 == "P", c(1, 3)], all.x = T)
kdata2$t0 <- as.numeric(as.Date(as.character(kdata2$postvax), format = "%d/%m/%Y") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1
kdata2$T <- as.numeric(as.Date(as.character(kdata2$date), format = "%Y%m%d") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1
kdata2$T[is.na(kdata2$T)] <- as.numeric(as.Date(as.character(kdata2$end.date[is.na(kdata2$T)]), format = "%d/%m/%Y") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1

kdata3 <- kdata2[kdata2$t0 < kdata2$T, ]

w1 <- model.p1$coef[2] # waning parameter run waning script first
w0 <- model.p0$coef[2]

kdata3.p1 <- kdata3[kdata3$intervention == "TIV", ]
kdata3.p0 <- kdata3[kdata3$intervention == "placebo", ]

mrow <- kdata3.p1$T[1] - kdata3.p1$t0[1] + 1
kdata.td.p1 <- data.frame(
  hhID = rep(kdata3.p1$hhID[1], mrow), TIV = 1 * (kdata3.p1$intervention[1] == "TIV"), postvax.P = log2(kdata3.p1$postvax.pH1[1]) + w1 * ((1:mrow) - 1), t = (1:mrow) - 1,
  event = c(rep(0, mrow - 1), kdata3.p1$swab.pH1[1])
)
kdata.td.p1$proxy <- s.proxy$s.proxy[kdata3.p1$t0[1]:kdata3.p1$T[1]]
kdata.td2.p1 <- kdata.td.p1[nrow(kdata.td.p1), ]

for (j in 2:nrow(kdata3.p1)) {
  mrow <- kdata3.p1$T[j] - kdata3.p1$t0[j] + 1
  temp <- data.frame(
    hhID = rep(kdata3.p1$hhID[j], mrow), TIV = 1 * (kdata3.p1$intervention[j] == "TIV"), postvax.P = log2(kdata3.p1$postvax.pH1[j]) + w1 * ((1:mrow) - 1), t = (1:mrow) - 1,
    event = c(rep(0, mrow - 1), kdata3.p1$swab.pH1[j])
  )
  temp$proxy <- s.proxy$s.proxy[kdata3.p1$t0[j]:kdata3.p1$T[j]]
  temp2 <- temp[nrow(temp), ]
  kdata.td.p1 <- rbind(kdata.td.p1, temp)
  kdata.td2.p1 <- rbind(kdata.td2.p1, temp2)
}

mrow <- kdata3.p0$T[1] - kdata3.p0$t0[1] + 1
kdata.td.p0 <- data.frame(
  hhID = rep(kdata3.p0$hhID[1], mrow), TIV = 1 * (kdata3.p0$intervention[1] == "TIV"), postvax.P = log2(kdata3.p0$postvax.pH1[1]) + w0 * ((1:mrow) - 1), t = (1:mrow) - 1,
  event = c(rep(0, mrow - 1), kdata3.p0$swab.pH1[1])
)
kdata.td.p0$proxy <- s.proxy$s.proxy[kdata3.p0$t0[1]:kdata3.p0$T[1]]
kdata.td2.p0 <- kdata.td.p0[nrow(kdata.td.p0), ]

for (j in 2:nrow(kdata3.p0)) {
  mrow <- kdata3.p0$T[j] - kdata3.p0$t0[j] + 1
  temp <- data.frame(
    hhID = rep(kdata3.p0$hhID[j], mrow), TIV = 1 * (kdata3.p0$intervention[j] == "TIV"), postvax.P = log2(kdata3.p0$postvax.pH1[j]) + w0 * ((1:mrow) - 1), t = (1:mrow) - 1,
    event = c(rep(0, mrow - 1), kdata3.p0$swab.pH1[j])
  )
  temp$proxy <- s.proxy$s.proxy[kdata3.p0$t0[j]:kdata3.p0$T[j]]
  temp2 <- temp[nrow(temp), ]
  kdata.td.p0 <- rbind(kdata.td.p0, temp)
  kdata.td2.p0 <- rbind(kdata.td2.p0, temp2)
}

kdata.td.p <- rbind(kdata.td.p0, kdata.td.p1)

# This model fit is where the figure derives its coefficients
cox.td.A.m <- coxph(Surv(t, event) ~ postvax.P + proxy, data = kdata.td.p)

#---------------------------------------------------------------------------------------------------------------------------------

#
# Influenza B  Model
#

kdata <- read.csv(paste(dir, "serology.csv", sep = ""))
demog <- read.csv(paste(dir, "demographic.csv", sep = ""))
swab <- read.csv(paste(dir, "swab.csv", sep = ""))
random <- read.csv(paste(dir, "randomcode.csv", sep = ""))
fluproxy <- read.csv("http://web.hku.hk/~bcowling/influenza/data/influenza_proxy_1998to2013.csv")

fluproxy <- fluproxy[614:679, ]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending, format = "%d/%m/%Y") - as.Date(time0, format = "%d/%m/%Y")) - 1
proxy <- data.frame(ti = fluproxy$ti, proxy = (1 / r) * fluproxy$B.proxy)

# create smoothed proxy for background risk
s.proxy <- data.frame(ti = 1:max(proxy$ti), s.proxy = predict(smooth.spline(proxy$ti, proxy$proxy, spar = 0.5), 1:max(proxy$ti))$y)

kdata <- kdata[kdata$member == 0, ]
kdata <- merge(kdata, random, by = "hhID", all.x = T)
kdata <- merge(kdata, demog[c(1, 2, 15)], by = c("hhID", "member"), all.x = T)

swab <- swab[swab$FluB == "P", ]
swab <- swab[order(swab$hhID, swab$member, swab$date), ]
swab$mark <- 0
for (i in 2:nrow(swab)) {
  if (swab$hhID[i] == swab$hhID[i - 1] & swab$member[i] == swab$member[i - 1]) swab$mark[i] <- 1
}
swab <- swab[swab$mark == 0 & swab$member == 0, ]

# merge in PCR type and date
kdata$swab.Victoria <- 1 * (kdata$hhID %in% swab$hhID[swab$FluB.subtype == "Victoria"])
kdata$swab.Yamagata <- 1 * (kdata$hhID %in% swab$hhID[swab$FluB.subtype == "Yamagata"])
kdata <- kdata[!is.na(kdata$postvax.B.Brisbane) & kdata$postvax != "", ]
kdata2 <- merge(kdata, swab[swab$FluB.subtype == "Victoria" | swab$FluB.subtype == "Yamagata", c(1, 3)], all.x = T)
kdata2$t0 <- as.numeric(as.Date(as.character(kdata2$postvax), format = "%d/%m/%Y") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1
kdata2$T <- as.numeric(as.Date(as.character(kdata2$date), format = "%Y%m%d") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1
kdata2$T[is.na(kdata2$T)] <- as.numeric(as.Date(as.character(kdata2$end.date[is.na(kdata2$T)]), format = "%d/%m/%Y") - as.Date(as.character("1/10/2009"), format = "%d/%m/%Y")) + 1

kdata3 <- kdata2[kdata2$t0 < kdata2$T, ]

# construct a dataframe for best waning paramter value

w1 <- c(model.b1.age$coef[2], model.b1.age$coef[4]) # waning parameter run waning script first
w0 <- model.b0$coef[2]

kdata3.b1 <- kdata3[kdata3$intervention == "TIV", ]
kdata3.b0 <- kdata3[kdata3$intervention == "placebo", ]

mrow <- kdata3.b1$T[1] - kdata3.b1$t0[1] + 1
kdata.td.b1 <- data.frame(
  hhID = rep(kdata3.b1$hhID[1], mrow), TIV = 1 * (kdata3.b1$intervention[1] == 1), postvax.b = log2(kdata3.b1$postvax.B.Brisbane[1]) + w1[1] * ((1:mrow) - 1) + w1[2] * 1 * (kdata3.b1$age[1] < 9) * ((1:mrow) - 1), t = (1:mrow) - 1,
  event = c(rep(0, mrow - 1), 1 * (kdata3.b1$swab.Victoria[1] == 1))
)
kdata.td.b1$proxy <- s.proxy$s.proxy[kdata3.b1$t0[1]:kdata3.b1$T[1]]
kdata.td2.b1 <- kdata.td.b1[nrow(kdata.td.b1), ]

for (j in 2:nrow(kdata3.b1)) {
  mrow <- kdata3.b1$T[j] - kdata3.b1$t0[j] + 1
  temp <- data.frame(
    hhID = rep(kdata3.b1$hhID[j], mrow), TIV = 1 * (kdata3.b1$intervention[j] == 1), postvax.b = log2(kdata3.b1$postvax.B.Brisbane[j]) + w1[1] * ((1:mrow) - 1) + w1[2] * 1 * (kdata3.b1$age[j] < 9) * ((1:mrow) - 1), t = (1:mrow) - 1,
    event = c(rep(0, mrow - 1), 1 * (kdata3.b1$swab.Victoria[j] == 1))
  )
  temp$proxy <- s.proxy$s.proxy[kdata3.b1$t0[j]:kdata3.b1$T[j]]
  temp2 <- temp[nrow(temp), ]
  kdata.td.b1 <- rbind(kdata.td.b1, temp)
  kdata.td2.b1 <- rbind(kdata.td2.b1, temp2)
}

mrow <- kdata3.b0$T[1] - kdata3.b0$t0[1] + 1
kdata.td.b0 <- data.frame(
  hhID = rep(kdata3.b0$hhID[1], mrow), TIV = 1 * (kdata3.b0$intervention[1] == 1), postvax.b = log2(kdata3.b0$postvax.B.Brisbane[1]) + w0 * ((1:mrow) - 1), t = (1:mrow) - 1,
  event = c(rep(0, mrow - 1), 1 * (kdata3.b0$swab.Victoria[1] == 1))
)
kdata.td.b0$proxy <- s.proxy$s.proxy[kdata3.b0$t0[1]:kdata3.b0$T[1]]
kdata.td2.b0 <- kdata.td.b0[nrow(kdata.td.b0), ]

for (j in 2:nrow(kdata3.b0)) {
  mrow <- kdata3.b0$T[j] - kdata3.b0$t0[j] + 1
  temp <- data.frame(
    hhID = rep(kdata3.b0$hhID[j], mrow), TIV = 1 * (kdata3.b0$intervention[j] == 1), postvax.b = log2(kdata3.b0$postvax.B.Brisbane[j]) + w0 * ((1:mrow) - 1), t = (1:mrow) - 1,
    event = c(rep(0, mrow - 1), 1 * (kdata3.b0$swab.Victoria[j] == 1))
  )
  
  temp$proxy <- s.proxy$s.proxy[kdata3.b0$t0[j]:kdata3.b0$T[j]]
  temp2 <- temp[nrow(temp), ]
  kdata.td.b0 <- rbind(kdata.td.b0, temp)
  kdata.td2.b0 <- rbind(kdata.td2.b0, temp2)
}

kdata.td.b <- rbind(kdata.td.b0, kdata.td.b1)

vic.vic <- cox.td.B <- coxph(Surv(t, event) ~ postvax.b + proxy, data = kdata.td.b)

# My script ===================================================================

all_soph_dat <- bind_rows(
  mod_soph_dat(kdata.td.p, "A H1N1pdm09"), 
  mod_soph_dat(kdata.td.b, "B Vic")
)

save_data(all_soph_dat, "sophia")

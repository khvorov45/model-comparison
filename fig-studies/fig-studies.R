# Creates figure showing study designs
# Sheena Sullivan

# Directories used
fig_studies_dir <- here::here("fig-studies")

# Script ======================================================================

pdf(file = file.path(fig_studies_dir, "fig-studies.pdf"), 8, 3.5)
par(mai = c(0, 2, 0, 0))
plot(
  NA, NA,
  xlim = c(0, 5), ylim = c(0, 4.5), axes = F,
  xlab = "Study timeline", ylab = ""
)
axis(
  2,
  at = seq(0.5, 3.5, 1),
  labels = rep(c("Observational", "Challenge"), 2),
  las = 2,
  tick = F,
  cex.axis = 0.8
)
mtext(
  "Vaccinated \nsubjects",
  side = 2,
  cex = 1,
  lwd = 2,
  font = 2,
  line = 7,
  at = c(3, -5)
)
mtext(
  "Unvaccinated \nsubjects",
  side = 2,
  cex = 1,
  lwd = 2,
  font = 2,
  line = 7,
  at = c(1, -5)
)
legend(
  "topright",
  pch = c(16, 6, 13, 14), col = c(2, 6, 1, 3), box.lty = 0,
  c("Blood draw", "Vaccination", "Challenge", "Swab"), horiz = T, cex = 0.7
)
legend(
  2.4, 4.4,
  lty = c(1, 1), col = c(4, 3), box.lty = 0,
  c("Challenge infection", "Seasonal epidemic"), horiz = T, cex = 0.7
)

# pre-vaccination bleed
points(0, 2.12, pch = 16, col = 2, cex = 1.2)
points(0, 3.12, pch = 16, col = 2, cex = 1.2)
# vaccinations
vpos <- 0.075
text(vpos, c(2.5, 3.5), "Vaccination", cex = 0.5)
arrows(vpos, 2.45, vpos, 2.2, length = 0.05)
arrows(vpos, 3.45, vpos, 3.2, length = 0.05)
points(vpos, 2.13, pch = 6, cex = 1, col = 6)
points(vpos, 3.13, pch = 6, cex = 1, col = 6)

# post-vaccination or pre-season/challenge bleed
for (i in 0:3) {
  points(0.5, .12 + i, pch = 16, col = 2, cex = 1.2)
}
# pre-season bleed
pspos <- 0.5
text(pspos, c(0.6, 2.6), "Pre-\nseason", cex = 0.5)
arrows(pspos, 0.5, pspos, .2, length = 0.05)
arrows(pspos, 2.5, pspos, 2.2, length = 0.05)

# challenge
cpos <- 0.575
text(cpos, c(1.5, 3.5), "Challenge", cex = 0.5)
arrows(cpos, 1.45, cpos, 1.2, length = 0.05)
arrows(cpos, 3.45, cpos, 3.2, length = 0.05)
points(cpos, 1.12, pch = 13, cex = 1.3, col = 1)
points(cpos, 3.12, pch = 13, cex = 1.3, col = 1)
# challenge epidemic:
x <- seq(0, 5, length = 1000)
y <- dlnorm(x, mean = 0.02, sd = 0.08)
lines(x - 0.15, y / 20 + 1, lwd = 1, col = 4)
lines(x - 0.15, y / 20 + 3, lwd = 1, col = 4)

# challenge swab
cspos <- 0.85
text(cspos, c(1.75, 3.75), "Swab", cex = 0.5)
arrows(cspos, 1.65, cspos, 1.35, length = 0.05)
arrows(cspos, 3.65, cspos, 3.35, length = 0.05)
points(cspos, 1.12, pch = 14, cex = 1.3, col = 3)
points(cspos, 3.12, pch = 14, cex = 1.3, col = 3)

# post challenge bleed
text(1.1, c(1.6, 3.6), "Post \nChallenge", cex = 0.5)
arrows(1.1, 1.45, 1.1, 1.2, length = 0.05)
arrows(1.1, 3.45, 1.1, 3.2, length = 0.05)
points(1.1, 1.1, pch = 16, col = 2, cex = 1.2)
points(1.1, 3.1, pch = 16, col = 2, cex = 1.2)

# natural epidemic
curve(dnorm(x, mean = 2.5, sd = 0.5) / 2, add = TRUE, col = 3)
curve(dnorm(x, mean = 2.5, sd = 0.5) / 2 + 1.9999, add = TRUE, col = 3)
# post-season bleed
pspos <- 4
text(pspos, c(0.6, 2.6), "Post-\nseason", cex = 0.5)
arrows(pspos, 0.45, pspos, .2, length = 0.05)
arrows(pspos, 2.45, pspos, 2.2, length = 0.05)
points(4, 0.12, pch = 16, col = 2, cex = 1.2)
points(4, 2.12, pch = 16, col = 2, cex = 1.2)

# lines between designs
abline(h = c(0:4))

# outcomes
text(4.5, c(0.5, 1.5), "PCR or \nseroconversion", cex = 0.6)
text(4.5, c(2.5, 3.5), "PCR", cex = 0.6)

dev.off()

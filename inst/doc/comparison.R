## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.5)

## ----setup--------------------------------------------------------------------
library(pcdreg)

## ----bothmodels---------------------------------------------------------------
set.seed(1)
d <- sim_pcd(300, beta = c(1, -1), lambda = function(t) 8 / (1 + t))

rate <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
mean <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
               model = "mean")

round(cbind(truth = c(1, -1), rate = coef(rate), mean = coef(mean)), 3)

## ----monotone-----------------------------------------------------------------
mu <- baseline(mean)$mean
c(increasing = !is.unsorted(mu), decreases_at = sum(diff(mu) < 0))

## ----predmono-----------------------------------------------------------------
pr <- predict(rate, d, type = "mean")
all(vapply(split(pr$mean, pr$id), function(v) all(diff(v) >= 0), logical(1)))

## ----layout-------------------------------------------------------------------
head(subset(d, id == 1), 4)

## ----robust-------------------------------------------------------------------
set.seed(2)
od <- sim_pcd(300, beta = c(1, -1), frailty = 1)  # gamma frailty
odfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = od)

round(rbind(robust      = sqrt(diag(vcov(odfit, "robust"))),
            information = sqrt(diag(vcov(odfit, "information")))), 4)

## ----dataplot, fig.alt = "Panel count data drawn as tiles: one row per subject, one tile per examination interval spanning the time it covers, shaded from pale to dark with the number of events recorded in that interval."----
autoplot(with(d, pcd(id, tstart, tstop, count)), max_subjects = 40)

## ----bench-code, eval = FALSE-------------------------------------------------
# sim_fixed <- function(n, beta = c(1, -1), tau = 2, exam_mean = 4) {
#   Lambda <- function(t) 8 * log(1 + t)
#   do.call(rbind, lapply(seq_len(n), function(i) {
#     nexam <- max(1L, rpois(1, exam_mean))
#     exams <- sort(unique(round(runif(nexam, 0, runif(1, 0.9 * tau, tau)), 2)))
#     exams <- exams[exams > 0]
#     x1 <- runif(1); x2 <- runif(1)
#     mu <- diff(c(0, Lambda(exams))) * exp(beta[1] * x1 + beta[2] * x2)
#     data.frame(id = i, time = exams, count = rpois(length(exams), mu),
#                x1 = x1, x2 = x2)
#   }))
# }
# 
# d <- sim_fixed(200)
# 
# # pcdreg: the standard error comes back with the fit.
# system.time({
#   f <- pcdreg(pcd(id, time, count) ~ x1 + x2, data = d, model = "mean")
#   sqrt(diag(vcov(f, "robust")))
# })
# 
# # spef: a point estimate, then a standard error by resampling.
# system.time(panelReg(PanelSurv(id, time, count) ~ x1 + x2, d,
#                      method = "EE.HSWc", se = "NULL"))
# system.time(panelReg(PanelSurv(id, time, count) ~ x1 + x2, d,
#                      method = "EE.HSWc", se = "Bootstrap"))

## ----bench-data---------------------------------------------------------------
bench <- data.frame(
  method = rep(c("pcdreg mean", "pcdreg rate", "spef AEE", "spef EE.HSWc",
                 "spef EE.HSWm", "spef MPL"), each = 4, times = 2),
  framing = rep(c("estimate + SE", "estimate only"), each = 24),
  n = rep(c(50, 100, 200, 400), times = 12),
  seconds = c(
    0.00, 0.00, 0.00, 0.00,   0.66, 1.83, 4.84, 7.23,
    0.89, 1.65, 3.17, 6.38,   0.15, 0.35, 0.79, 3.53,
    0.45, 0.89, 1.77, 5.23,   1.57, 2.27, 3.69, 4.98,
    0.00, 0.00, 0.00, 0.00,   0.66, 1.80, 4.89, 7.30,
    0.08, 0.11, 0.08, 0.11,   0.01, 0.03, 0.04, 0.14,
    0.01, 0.05, 0.08, 0.19,   0.04, 0.11, 0.12, 0.24
  )
)
knitr::kable(
  stats::reshape(bench, idvar = c("method", "framing"), timevar = "n",
                 direction = "wide"),
  col.names = c("method", "framing", "n = 50", "n = 100", "n = 200", "n = 400"),
  row.names = FALSE
)

## ----bench-plot, fig.alt = "Two panels of median run time against sample size. With standard errors included, pcdreg's mean model is flat at essentially zero while every spef method rises to between 3.5 and 6.4 seconds at n = 400. For point estimates alone the ordering reverses: spef stays under a quarter of a second and pcdreg's rate model rises to 7.3 seconds.", fig.height = 3.6----
library(ggplot2)
ggplot(bench, aes(n, seconds, colour = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  facet_wrap(~ framing) +
  scale_colour_manual(values = c(
    "pcdreg mean"  = "#3987e5", "pcdreg rate"  = "#0d366b",
    "spef AEE"     = "#b0aca4", "spef EE.HSWc" = "#8a857c",
    "spef EE.HSWm" = "#645f57", "spef MPL"     = "#3f3b35")) +
  labs(x = "Subjects", y = "Median seconds", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())


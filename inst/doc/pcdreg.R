## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(pcdreg)

## ----three-arg----------------------------------------------------------------
visits <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 0, 3, 1, 4),
  dose  = c(10, 10, 10, 25, 25)
)
pcd(visits$id, visits$time, visits$count)

## ----four-arg-----------------------------------------------------------------
cp <- data.frame(
  id     = c(1, 1, 1, 2),
  tstart = c(0.0, 0.5, 0.9, 0.0),
  tstop  = c(0.5, 0.9, 1.2, 0.8),
  count  = c(2, NA, 0, 3),
  dose   = c(10, 10, 25, 25)
)
pcd(cp$id, cp$tstart, cp$tstop, cp$count)

## ----convert, eval = FALSE----------------------------------------------------
# new <- old
# new$count[old$event == 0] <- NA
# pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = new)

## ----fit----------------------------------------------------------------------
set.seed(1)
d <- sim_pcd(200, beta = c(1, -1), lambda = function(t) 8 / (1 + t))
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
summary(fit)

## ----variance-----------------------------------------------------------------
se <- function(type) sqrt(diag(vcov(fit, type)))
rbind(robust = se("robust"), information = se("information"))

## ----matrices-----------------------------------------------------------------
fit$Omega
fit$S

## ----profile------------------------------------------------------------------
pfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d,
                  profile = TRUE)
rbind(information = sqrt(diag(vcov(pfit, "information"))),
      profile     = sqrt(diag(vcov(pfit, "profile"))))

## ----overdispersed------------------------------------------------------------
set.seed(2)
od <- sim_pcd(200, beta = c(1, -1), frailty = 1)
odfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = od)
rbind(robust = sqrt(diag(vcov(odfit, "robust"))),
      information = sqrt(diag(vcov(odfit, "information"))))

## ----baseline, fig.alt = "Estimated cumulative baseline rate rising with time, with the true curve overlaid as a dashed red line and closely tracking it."----
head(baseline(fit))

truth <- data.frame(time = fit$baseline$time,
                    value = 8 * log(1 + fit$baseline$time))
autoplot(fit) +
  ggplot2::geom_line(data = truth, colour = "#c0392b", linetype = 2)

## ----predict, fig.alt = "Predicted mean count for one subject as a step function, rising monotonically over follow-up."----
pred <- predict(fit, d, type = "mean")
head(pred)

one <- pred[pred$id == 1, ]
ggplot2::ggplot(one, ggplot2::aes(x = time, y = mean)) +
  ggplot2::geom_step(colour = "#2a78d6", linewidth = 0.6) +
  ggplot2::labs(x = "Time", y = "Predicted mean count", title = "Subject 1") +
  ggplot2::theme_minimal(base_size = 11)

## ----tidiers------------------------------------------------------------------
tidy(fit, conf.int = TRUE)
glance(fit)

## ----exponentiate-------------------------------------------------------------
tidy(fit, conf.int = TRUE, exponentiate = TRUE)

## ----augment------------------------------------------------------------------
aug <- augment(fit)
head(aug[!is.na(aug$.observed),
         c("id", "tstop", ".fitted", ".observed", ".resid")])

## ----meanmodel----------------------------------------------------------------
mfit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d, model = "mean")
summary(mfit)

## ----compare------------------------------------------------------------------
cbind(rate = coef(fit), mean = coef(mfit))

## ----monotone-----------------------------------------------------------------
mu <- baseline(mfit)$mean
c(increasing = !is.unsorted(mu), decreases_at = sum(diff(mu) < 0))

## ----passes-------------------------------------------------------------------
c(iterations = fit$iterations, EM_passes = fit$passes,
  converged = fit$converged)


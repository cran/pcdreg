## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(pcdreg)

## ----three-arg----------------------------------------------------------------
visits <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 0, 3, 1, 4),
  age   = c(24, 24, 24, 31, 31),
  arm   = c("drug", "drug", "drug", "placebo", "placebo")
)

pcd(visits$id, visits$time, visits$count)

## ----three-arg-fit, eval = FALSE----------------------------------------------
# pcdreg(pcd(id, time, count) ~ age + arm, data = visits)

## ----four-arg-----------------------------------------------------------------
cp <- data.frame(
  id     = c(1, 1, 1, 2),
  tstart = c(0.0, 0.5, 0.9, 0.0),
  tstop  = c(0.5, 0.9, 1.2, 0.8),
  count  = c(2, NA, 0, 3),
  dose   = c(10, 10, 25, 25)
)

pcd(cp$id, cp$tstart, cp$tstop, cp$count)

## ----rhs, eval = FALSE--------------------------------------------------------
# pcd(id, tstart, tstop, count) ~ x1 + x2          # main effects
# pcd(id, tstart, tstop, count) ~ x1 * x2          # interaction
# pcd(id, tstart, tstop, count) ~ x1 + log(x2)     # transformation
# pcd(id, tstart, tstop, count) ~ x1 + factor(site)
# pcd(id, tstart, tstop, count) ~ .                # everything else

## ----factors------------------------------------------------------------------
set.seed(1)
d <- sim_pcd(60)
d$site <- factor(c("north", "south", "east")[1 + d$id %% 3])

fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + site, data = d)
names(coef(fit))

## ----dot----------------------------------------------------------------------
names(coef(pcdreg(pcd(id, tstart, tstop, count) ~ ., data = d)))

## ----cumulative---------------------------------------------------------------
# A frame whose counts accumulate: 2, 2, 5 rather than 2, 0, 3.
cumdat <- data.frame(
  id    = c(1, 1, 1, 2, 2),
  time  = c(0.5, 1.1, 1.9, 0.7, 1.4),
  count = c(2, 2, 5, 1, 5),
  age   = c(24, 24, 24, 31, 31)
)

cumdat$increment <- ave(cumdat$count, cumdat$id,
                        FUN = function(v) v - c(0, v[-length(v)]))
cumdat[, c("id", "time", "count", "increment")]

## ----cumulative-cp, eval = FALSE----------------------------------------------
# at <- !is.na(dat$count)
# dat$count[at] <- ave(dat$count[at], dat$id[at],
#                      FUN = function(v) v - c(0, v[-length(v)]))

## ----legacy, eval = FALSE-----------------------------------------------------
# new <- old
# new$count[old$event == 0] <- NA
# pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = new)

## ----merge--------------------------------------------------------------------
exams <- data.frame(id = c(1, 1), time = c(0.5, 1.2), count = c(2, 0))
covar <- data.frame(id = 1, time = 0.9, dose = 25)
start <- data.frame(id = 1, time = 0, dose = 10)

both <- merge(
  data.frame(id = c(exams$id, covar$id), tstop = c(exams$time, covar$time),
             count = c(exams$count, NA)),
  rbind(start[, c("id", "time", "dose")], covar[, c("id", "time", "dose")]),
  by.x = c("id", "tstop"), by.y = c("id", "time"), all.x = TRUE
)
both <- both[order(both$id, both$tstop), ]
both$dose <- ave(both$dose, both$id,
                 FUN = function(v) { for (i in seq_along(v))
                   if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]; v })
both$tstart <- ave(both$tstop, both$id,
                   FUN = function(v) c(0, v[-length(v)]))
both[, c("id", "tstart", "tstop", "count", "dose")]

## ----err-origin, error = TRUE-------------------------------------------------
try({
bad <- data.frame(id = 1, tstart = 0.2, tstop = 1, count = 3, x = 1)
pcdreg(pcd(id, tstart, tstop, count) ~ x, data = bad)
})

## ----err-gap, error = TRUE----------------------------------------------------
try({
bad <- data.frame(id = c(1, 1), tstart = c(0, 0.7), tstop = c(0.5, 1),
                  count = c(2, 1), x = c(1, 1))
pcdreg(pcd(id, tstart, tstop, count) ~ x, data = bad)
})

## ----err-count, error = TRUE--------------------------------------------------
try({
pcd(1:2, c(1, 2), c(0.5, 1.2))
})

## ----dropped------------------------------------------------------------------
# A censoring row appended after subject 1's final visit.
last1 <- max(d$tstop[d$id == 1])
trailing <- transform(d[d$id == 1, ][1, ],
                      tstart = last1, tstop = last1 + 0.5, count = NA)

fit2 <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2,
                  data = rbind(d, trailing))

## ----grid---------------------------------------------------------------------
fit <- pcdreg(pcd(id, tstart, tstop, count) ~ x1 + x2, data = d)
c(subjects = fit$n, examinations = fit$nexam,
  grid_times = fit$ngrid, events = fit$nevent)

## ----sanity-------------------------------------------------------------------
c(rows_supplied = nrow(d),
  examination_rows = sum(!is.na(d$count)),
  events = sum(d$count, na.rm = TRUE))


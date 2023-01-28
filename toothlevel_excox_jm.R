library(dplyr)
library(survival)
library(Matrix)
library(lme4)
library(nlme)
library(splines)
library(statmod)
library(JointModel)
library(MASS)
library(JM)
library(Rcpp)
library(here)
library(tidyverse)
library(kableExtra)
library(ggridges)
library(rstanarm)

# read the toothlevel data
vadls <- read.csv("vadls_jm_toothlevel.csv", header = TRUE)

## longitudinal data
valong <- vadls

## survival data
vasurv <- vadls %>%
  dplyr::group_by(id, tooth) %>%
  mutate(basepock = maxpock[1L],
         basesmoke = smoking[1L]) %>%
  arrange(obstime) %>%
  slice(n()) %>%
  ungroup() 

## data for extended Cox model 
tdsurv <- tmerge(vasurv, vasurv, id = idtooth, endpt = event(time, tloss) )
tdsurv <- tmerge(tdsurv, vadls, id = idtooth,
                 maxpock = tdc(obstime, maxpock),
                 smoking = tdc(obstime, smoking),
                 bmi = tdc(obstime, basebmi))


# extended Cox model
tdcox <- coxph(Surv(tstart, tstop, endpt) ~ maxpock + baseage + basebmi + basesmoke +college,
               data = tdsurv, cluster = id)

summary(tdcox)



# joint model
# first, fit longitudinal submodel
jmlong <- lme(maxpock ~ obstime,
              random = ~ obstime | idtooth,
              data = valong,
              na.action = na.omit)
summary(jmlong)

# next, fit survival submodel
jmcox <- coxph(Surv(time, tloss) ~ baseage + basebmi + basesmoke +college,
               cluster = id,
               data = vasurv,
               x = TRUE)
summary(jmcox)

# finally, fit joint model (current value model)
jmtloss <- jointModel(jmlong, jmcox,
                      timeVar = "obstime",
                      method = "piecewise-PH-aGH")

summary(jmtloss)


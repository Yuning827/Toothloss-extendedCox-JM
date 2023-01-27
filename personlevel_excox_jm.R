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

# read the personlevel data
tldat <- read.csv("vadls_jm_personlevel.csv", header = TRUE)

## data for extended Cox model 
tddat <- tldat %>%
  rename(tstart = year) %>%
  group_by(id) %>%
  mutate(tstop = dplyr::lead(tstart),
         TL = dplyr::lead(TL)) %>%
  dplyr::filter(row_number()!=n()|n()==1)

## data for joint model
tlLong <- tldat 

lastdat <- tldat %>%
  group_by(id) %>%
  arrange(year) %>%
  dplyr::slice(n()) %>%
  ungroup 

tlSurv <- lastdat %>%
  dplyr::select(id, year, TL, baseage, college, basesmoke, basebmi, basenumteeth) %>%
  rename(years=year)  


### extended Cox model
td.tl <- coxph(Surv(tstart, tstop, TL) ~ pctpocket5mm + pctabl40 + pctmobil05mm + 
                 baseage + basebmi + basesmoke +college + basenumteeth, data = tddat)
summary(td.tl)  


### joint model

jm.tl <- stan_jm(formulaLong = list(pctpocket5mm ~ year + (1 + year| id), # adding baseline covariates did not change the results
                                    pctmobil05mm ~ year + (1 + year| id),
                                    pctabl40 ~ year + (1 + year| id)),
                 dataLong = tlLong,
                 formulaEvent = Surv(years, TL) ~ baseage + basebmi + basesmoke +college + basenumteeth,
                 dataEvent = tlSurv,
                 time_var = "year",
                 # this next line is only to keep the example small in size!
                 chains = 1, cores = 1, seed = 12345, iter = 1000)
print(jm.tl)
summary(jm.tl, probs = c(0.025, 0.975))
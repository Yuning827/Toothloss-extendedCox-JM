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
#library(webshot)
library(kableExtra)
library(ggridges)
library(rstanarm)

tldat <- read.csv("vadls_jm_person_synthetic.csv", header = TRUE)



lastdat <- tldat %>%
  group_by(id) %>%
  arrange(year) %>%
  dplyr::slice(n()) %>%
  ungroup 

### data for extended Cox model 
tddat <- tldat %>%
  rename(tstart = year) %>%
  group_by(id) %>%
  mutate(tstop = dplyr::lead(tstart),
         TL = dplyr::lead(TL)) %>%
  dplyr::filter(row_number()!=n()|n()==1)


### data for joint model
tlLong <- tldat %>%
  dplyr::select(id, year,  pctpocket5mm, pctmobil05mm, pctabl40, TL,
                baseage, college, basesmoke, basebmi, basenumteeth)
tlSurv <- lastdat %>%
  dplyr::select(id, year, TL, baseage, college, basesmoke, basebmi, basenumteeth) %>%
  rename(years=year)  

tldat2 <- tldat %>%
  dplyr::select(id, year,  pctpocket5mm, pctmobil05mm, pctabl40, TL,
                baseage, college, basesmoke, basebmi, basenumteeth)


length(unique(tlLong$id))
length(unique(tlSurv$id))
length(unique(tldat$id))
sum(tlSurv$basenumteeth)
